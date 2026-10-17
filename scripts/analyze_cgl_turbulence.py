"""Snapshot PDFs and spectra for uniform AthenaK grids, including uniform AMR leaves.

One snapshot at a time; small NPZ caches separate calculation from plotting.
Run with --self-check, or MANIFEST (parallel shards use SLURM_PROCID/SLURM_NTASKS).
"""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy import fft

REPO = Path.home() / 'athenak-cgl'
sys.path.insert(0, str(REPO / 'vis/python'))
from bin_convert import read_binary

WORKERS = int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))
FIELDS = ('dens', 'velx', 'vely', 'velz', 'eint', 'p_perp', 'bcc1', 'bcc2', 'bcc3')
GRADIENTS = ('parallel_parallel', 'parallel_perp', 'perp_parallel', 'perp_perp')


def metadata(path):
    """Read only the small header; never load a field to discover its time."""
    with Path(path).open('rb') as stream:
        assert stream.readline().strip() == b'Athena binary output version=1.1'
        count = int(stream.readline().split(b'=')[-1])
        pre = dict(line.decode().strip().split('=', 1)
                   for line in (stream.readline() for _ in range(count-1)))
        pre = {key.strip(): value.strip() for key, value in pre.items()}
        nvars = int(stream.readline().split(b'=')[-1])
        names = stream.readline().decode().split()[1:]
        assert len(names) == nvars
        size = int(stream.readline().split(b'=')[-1])
        header = stream.read(size)
        assert len(header) == size
    params, section = {}, None
    for line in header.decode().splitlines():
        line = line.split('#')[0].strip()
        if line.startswith('<'):
            section = line.strip('<>')
        elif '=' in line:
            key, value = line.split('=', 1)
            params[section + '/' + key.strip()] = value.strip()
    return dict(time=float(pre['time']), cycle=int(pre['cycle']), params=params,
                header_sha256=hashlib.sha256(header).hexdigest(),
                nvars=nvars, varbytes=int(pre['size of variable']))


def assemble(raw):
    """Tile complete, same-level blocks; reject mixed AMR rather than FFT it incorrectly."""
    logical = raw['mb_logical']
    assert len(set(logical[:, 3])) == 1, 'Mixed-level AMR needs a separate resampling policy'
    assert np.all(logical[:, :3] >= 0)
    blocks = logical[:, :3].max(axis=0) + 1
    assert np.prod(blocks) == len(logical) == len(set(map(tuple, logical)))
    block = np.array([raw[f'nx{i}_out_mb'] for i in (1, 2, 3)])
    assert all(raw[f'nx{i}_out_mb'] == raw[f'nx{i}_mb'] for i in (1, 2, 3))
    shape = tuple((blocks*block)[::-1])
    lengths = tuple(raw[f'x{i}max']-raw[f'x{i}min'] for i in (1, 2, 3))
    assert all(length > 0 for length in lengths)
    fields = {name: np.empty(shape, dtype=np.float32) for name in FIELDS}
    for m, (lx, ly, lz, level) in enumerate(logical):
        for axis, logical_index in enumerate((lx, ly, lz), start=1):
            width = lengths[axis-1]/blocks[axis-1]
            lower = raw[f'x{axis}min'] + logical_index*width
            np.testing.assert_allclose(raw['mb_geometry'][m,2*(axis-1):2*axis],
                                       [lower,lower+width], rtol=1e-6, atol=1e-8)
        section = tuple(slice(int(start), int(start+width))
                        for start, width in zip((np.array([lz,ly,lx])*block[::-1]), block[::-1]))
        for name in FIELDS:
            fields[name][section] = raw['mb_data'][name][m]
    for name, field in fields.items():
        assert np.isfinite(field).all(), name
    return fields, lengths


def fourier_grid(shape, lengths):
    """Physical shells; k is |k| Lx/(2pi). Include real-FFT multiplicities."""
    nz, ny, nx = shape
    lx, ly, lz = lengths
    kx = fft.rfftfreq(nx, d=lx/nx).astype(np.float32)
    ky = fft.fftfreq(ny, d=ly/ny).astype(np.float32)
    kz = fft.fftfreq(nz, d=lz/nz).astype(np.float32)
    dk = lx/max(lengths)
    radius = lx*np.sqrt(kz[:,None,None]**2+ky[None,:,None]**2+kx[None,None,:]**2)
    shells = np.floor(radius/dk + 0.5).astype(np.int32)
    weights = np.full(nx//2+1, 2, dtype=np.float32)
    weights[0] = 1
    if nx % 2 == 0:
        weights[-1] = 1
        kx[-1] = 0
    if ny % 2 == 0:
        ky[ny//2] = 0
    if nz % 2 == 0:
        kz[nz//2] = 0
    # A real spectral derivative has zero even-grid Nyquist multiplier.
    K = [2*np.pi*kx[None,None,:], 2*np.pi*ky[None,:,None], 2*np.pi*kz[:,None,None]]
    bins = int(shells.max())+1
    return dict(shells=shells, weights=weights, K=K, k2=sum(k*k for k in K),
                dk=dk, bins=bins, k=np.arange(bins)*dk, shape=shape)


def bin_power(transformed, grid, factor=1):
    power = (transformed.real**2+transformed.imag**2)*factor
    power *= grid['weights'][None,None,:]
    return np.bincount(grid['shells'].ravel(), weights=power.ravel(),
                       minlength=grid['bins'])/grid['dk']


def spectrum(field, grid, centered=True):
    values = field-np.mean(field, dtype=np.float64) if centered else field
    transformed = fft.rfftn(values, norm='forward', workers=WORKERS)
    result = bin_power(transformed, grid)
    expected = np.mean(values.astype(np.float64)**2) if field.size < 100000 else np.mean(values*values, dtype=np.float64)
    assert np.isclose(result.sum()*grid['dk'], expected, rtol=4e-5, atol=1e-12)
    return result


def gradient_spectra(velocity, b, grid):
    """Project G_ij=partial_j u_i AFTER differentiation; do not differentiate b."""
    shape = velocity[0].shape
    gradient = np.empty((3,3)+shape, dtype=np.float32)
    total = np.zeros(grid['bins'])
    divergence = np.zeros(shape, dtype=np.float32)
    for i in range(3):
        transformed = fft.rfftn(velocity[i], norm='forward', workers=WORKERS)
        total += bin_power(transformed, grid, grid['k2'])
        for j in range(3):
            gradient[i,j] = fft.irfftn(1j*grid['K'][j]*transformed, s=shape,
                                       norm='forward', workers=WORKERS)
            if i == j:
                divergence += gradient[i,j]
    along = [sum(gradient[i,j]*b[j] for j in range(3)) for i in range(3)]
    component = [sum(b[i]*gradient[i,j] for i in range(3)) for j in range(3)]
    strain = sum(b[i]*along[i] for i in range(3))
    curves = dict(Pgrad=total, PS=spectrum(strain,grid), Pdiv=spectrum(divergence,grid),
                  PQ=spectrum(strain-divergence/3,grid))
    curves['gradient_total'] = 0.5*total
    curves['gradient_parallel_parallel'] = 0.5*spectrum(strain,grid,False)
    curves['gradient_parallel_perp'] = sum(0.5*spectrum(along[i]-b[i]*strain,grid,False) for i in range(3))
    curves['gradient_perp_parallel'] = sum(0.5*spectrum(component[j]-b[j]*strain,grid,False) for j in range(3))
    curves['gradient_perp_perp'] = sum(
        0.5*spectrum(gradient[i,j]-b[i]*component[j]-along[i]*b[j]+b[i]*b[j]*strain,grid,False)
        for i in range(3) for j in range(3))
    assert np.isclose(sum(curves['gradient_'+key].sum() for key in GRADIENTS),
                      curves['gradient_total'].sum(), rtol=4e-5, atol=1e-12)
    return curves


def analyze(path, nominal_time):
    start = time.monotonic()
    meta = metadata(path)
    estimate = Path(path).stat().st_size/(meta['nvars']*meta['varbytes'])*4*70
    assert estimate < 350*2**30, 'Whole-grid working-set estimate exceeds 350 GiB; use distributed FFT analysis'
    raw = read_binary(str(path))
    fields, lengths = assemble(raw)
    actual_time, cycle = float(raw['time']), int(raw['cycle'])
    del raw
    gc.collect()
    shape = fields['dens'].shape
    # Measured below on 768^3; this whole-domain CPU FFT route is not distributed.
    grid = fourier_grid(shape, lengths)
    rho, ppar, pperp = (fields[key] for key in ('dens','eint','p_perp'))
    velocity = [fields[key] for key in ('velx','vely','velz')]
    magnetic = [fields[key] for key in ('bcc1','bcc2','bcc3')]
    assert min(float(rho.min()),float(ppar.min()),float(pperp.min())) > 0
    bsqr = sum(b*b for b in magnetic)
    bmag = np.sqrt(bsqr)
    assert float(bmag.min()) > 0, 'Zero field has no local parallel direction'
    b0 = abs(float(meta['params']['problem/b0']))
    assert b0 > 0
    beta_delta = 2*(pperp-ppar)/bsqr
    b_edges = np.geomspace(1e-6,1e3,181)
    beta_edges = np.linspace(-3,3,241)
    b_hist = np.histogram(bmag/b0,b_edges)[0]
    beta_hist = np.histogram(beta_delta,beta_edges)[0]
    assert b_hist.sum() == bmag.size
    curves = dict(k=grid['k'], dk=grid['dk'], time=nominal_time, actual_time=actual_time,
                  shape=np.array(shape), lengths=np.array(lengths),
                  b_centers=np.sqrt(b_edges[:-1]*b_edges[1:]),
                  b_pdf=b_hist/(bmag.size*np.diff(b_edges)),
                  beta_centers=0.5*(beta_edges[:-1]+beta_edges[1:]),
                  beta_pdf=beta_hist/(bmag.size*np.diff(beta_edges)))
    stats = dict(time=nominal_time, actual_time=actual_time, cycle=cycle, shape=list(shape),
                 lengths=lengths, b0=b0, b_mean=float(np.mean(bmag,dtype=np.float64)),
                 b_std=float(np.std(bmag,dtype=np.float64)),
                 beta_pdf_outside_fraction=float(1-beta_hist.sum()/bmag.size),
                 beta_delta_min=float(beta_delta.min()), beta_delta_max=float(beta_delta.max()),
                 p_parallel_min=float(ppar.min()), p_perp_min=float(pperp.min()),
                 kinetic=float(0.5*np.mean(rho*sum(u*u for u in velocity),dtype=np.float64)),
                 magnetic=float(0.5*np.mean(bsqr,dtype=np.float64)),
                 u_rms=float(np.sqrt(np.mean(sum(u*u for u in velocity),dtype=np.float64))))
    del beta_delta
    curves['EK'] = sum(0.5*spectrum(np.sqrt(rho)*u,grid) for u in velocity)
    curves['EB'] = sum(0.5*spectrum(b,grid) for b in magnetic)
    for key, field in [('Pppar',ppar),('Ppperp',pperp),('PB2',bsqr/2)]:
        curves[key] = spectrum(field,grid)
    bhat = [b/bmag for b in magnetic]
    curves['Pupar'] = spectrum(sum(u*b for u,b in zip(velocity,bhat)),grid)
    curves.update(gradient_spectra(velocity,bhat,grid))
    band = (grid['k']>=4)&(grid['k']<=32)
    for key in ('EK','EB','Pppar','Ppperp','PB2','PS','PQ','Pdiv','Pgrad','Pupar',
                *('gradient_'+key for key in GRADIENTS)):
        stats[key+'_integral'] = float(curves[key].sum()*grid['dk'])
        stats[key+'_k4_32'] = float(curves[key][band].sum()*grid['dk'])
    # A pressure-balance correlation is more informative than amplitudes alone.
    p = pperp-np.mean(pperp,dtype=np.float64)
    m = bsqr/2-np.mean(bsqr/2,dtype=np.float64)
    covariance = np.mean(p*m,dtype=np.float64)
    denom = np.sqrt(np.mean(p*p,dtype=np.float64)*np.mean(m*m,dtype=np.float64))
    stats['pperp_magnetic_pressure_correlation'] = float(covariance/denom) if denom > 0 else None
    stats['elapsed_seconds'] = time.monotonic()-start
    stats['process_peak_rss_gib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20
    return curves, stats


def self_check():
    shape, lengths = (32,16,16),(1.,1.,2.)
    grid = fourier_grid(shape,lengths)
    z,y,x = np.meshgrid(np.arange(32)/16,np.arange(16)/16,np.arange(16)/16,indexing='ij')
    zero = np.zeros(shape,dtype=np.float32)
    wave = np.sin(np.pi*z).astype(np.float32)
    power = spectrum(wave,grid)
    assert np.isclose(power[1]*grid['dk'],0.5,rtol=3e-6)
    velocity = [(3*np.sin(6*np.pi*z)+4*np.sin(8*np.pi*y)).astype(np.float32),zero,
                (np.sin(2*np.pi*x)+2*np.sin(4*np.pi*z)).astype(np.float32)]
    curves = gradient_spectra(velocity,[zero,zero,zero+1],grid)
    for key,k,amplitude in zip(GRADIENTS,(2,3,1,4),(8,18,2,32)):
        expected = (amplitude*np.pi)**2/4
        assert np.isclose(curves['gradient_'+key][int(k/grid['dk'])]*grid['dk'],expected,rtol=4e-5),key
    b = [(0.6*np.cos(np.pi*z)).astype(np.float32),(0.6*np.sin(np.pi*z)).astype(np.float32),zero+0.8]
    varied = gradient_spectra(velocity,b,grid)
    assert not np.allclose(sum(varied['gradient_'+key] for key in GRADIENTS),varied['gradient_total'])
    # Refined logical extents must override the root-grid resolution.
    raw = dict(mb_logical=np.array([[x,0,0,3] for x in range(4)]),
               x1min=0.,x1max=1.,x2min=0.,x2max=1.,x3min=0.,x3max=1.)
    for i in (1,2,3): raw[f'nx{i}_mb']=raw[f'nx{i}_out_mb']=2
    raw['mb_geometry']=np.array([[m/4,(m+1)/4,0,1,0,1] for m in range(4)])
    raw['mb_data']={key:[np.full((2,2,2),m,dtype=np.float32) for m in range(4)] for key in FIELDS}
    fields,_=assemble(raw)
    assert fields['dens'].shape==(2,2,8) and fields['dens'][0,0,6]==3
    raw['mb_logical'][0,3]=2
    try: assemble(raw)
    except AssertionError: pass
    else: raise AssertionError('mixed-level grid accepted')
    # Odd-grid real FFT requires weight 2 on its last stored mode.
    odd=fourier_grid((9,9,9),(1.,1.,1.))
    q=np.broadcast_to(np.sin(8*np.pi*np.arange(9)/9),(9,9,9)).astype(np.float32)
    assert np.isclose(spectrum(q,odd).sum(),.5,rtol=3e-6)
    print('Rectangular modes, Parseval, four local projections, uniform AMR tiling, mixed-level rejection and odd FFT passed.',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest',nargs='?',type=Path)
    parser.add_argument('--self-check',action='store_true')
    parser.add_argument('--limit',type=int)
    args=parser.parse_args()
    if args.self_check:
        self_check()
        return
    manifest=json.loads(args.manifest.read_text())
    tasks=manifest['tasks']
    rank=int(os.environ.get('SLURM_PROCID','0'))
    ranks=int(os.environ.get('SLURM_NTASKS','1'))
    codehash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    tasks=tasks[rank::ranks]
    if args.limit is not None: tasks=tasks[:args.limit]
    for task in tasks:
        source=Path(task['source']); output=Path(task['cache']); output.parent.mkdir(parents=True,exist_ok=True)
        info=output.with_suffix('.json')
        source_stat=source.stat()
        provenance=dict(source=str(source),size=source_stat.st_size,mtime_ns=source_stat.st_mtime_ns,analyzer_sha256=codehash)
        if output.exists() and info.exists() and json.loads(info.read_text()).get('provenance')==provenance:
            print('Cached',output,flush=True)
            continue
        print('Starting',task,flush=True)
        curves,stats=analyze(source,task['time'])
        assert source.stat().st_size==source_stat.st_size and source.stat().st_mtime_ns==source_stat.st_mtime_ns,'Source changed during read'
        stats['provenance']=provenance
        temp=output.with_suffix('.tmp.npz')
        np.savez(temp,**curves); temp.replace(output)
        temp=info.with_suffix('.tmp.json'); temp.write_text(json.dumps(stats,indent=2,allow_nan=False)+'\n'); temp.replace(info)
        print('Completed',output,stats['elapsed_seconds'],'s peak',stats['process_peak_rss_gib'],'GiB',flush=True)
        del curves,stats
        gc.collect()


if __name__=='__main__': main()
