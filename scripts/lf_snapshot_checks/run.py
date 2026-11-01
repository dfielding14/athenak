"""Compare heat-flux closures at measured physical scales in one frozen active CGL snapshot."""
import argparse
import csv
import json
import shutil
import tempfile
from pathlib import Path
from time import perf_counter

import numpy as np

from fluxes import evaluate_faces
from optical import optical_faces
from scale_analysis import build_pairs, summarize
from snapshot_io import NativeSnapshot, load_fields


def stencils(pairs, shape):
    """Deduplicate the six faces per endpoint; lower faces are indexed by their right cell."""
    xyz, endpoint = np.unique(pairs.reshape(-1, 3), axis=0, return_inverse=True)
    faces = np.empty((len(xyz), 3, 2, 4), dtype=np.int64)
    for d in range(3):
        plus = xyz.copy()
        plus[:, d] = (plus[:, d] + 1) % shape[2-d]
        faces[:, d, :, 0] = d
        faces[:, d, 0, 1:] = xyz[:, ::-1]
        faces[:, d, 1, 1:] = plus[:, ::-1]
    ids, inv = np.unique(faces.reshape(-1, 4), axis=0, return_inverse=True)
    return xyz, endpoint.reshape(pairs.shape[:-1]), ids, inv.reshape(faces.shape[:-1])


def paired_values(fields, xyz, endpoint, inv, q, bn, bf, dx):
    samples = fields[xyz[:, 2], xyz[:, 1], xyz[:, 0]].astype(float)
    delta = (samples[:, 5] - samples[:, 4])[endpoint]
    B = np.linalg.norm(samples[:, :3], axis=1)
    if np.any(bf <= 0):
        raise ValueError('Undefined face magnetic moment')
    eflux = bn * (.5*q[..., 0] + q[..., 1])
    mflux = bn*q[..., 1]/bf
    # f[..., 0] is the low face and f[..., 1] is the high face.
    def minus_div(f):
        return ((f[:, inv[:, :, 0]] - f[:, inv[:, :, 1]])/dx).sum(axis=-1)
    tendency = (3*B*minus_div(mflux) - 2*minus_div(eflux))[:, endpoint]
    # Signed field-aligned scalars prior to projection; same six faces across methods.
    qcell = q[:, inv].mean(axis=(2, 3))[:, endpoint]
    return delta, tendency, qcell


def face_rows(fields, pairs, meta, result, groups, context, store, prefix):
    xyz, endpoint, ids, inv = stencils(pairs, fields.shape[:3])
    delta, tendency, qs = paired_values(fields, xyz, endpoint, inv, result['q'],
                                         result['bn'], result['bf'], np.array(context['spacing']))
    store[prefix+'_pairs_xyz'] = pairs
    store[prefix+'_face_ids_axis_zyx'] = ids
    for key in ('q','labels','bn','bf'):
        store[prefix+'_'+key] = np.asarray(result[key])
    store[prefix+'_paired_delta'] = delta
    store[prefix+'_paired_tendency'] = tendency
    displacement = pairs[...,1,:] - pairs[...,0,:]
    shape_xyz = np.array(fields.shape[:3][::-1])
    displacement = (displacement + shape_xyz//2) % shape_xyz - shape_xyz//2
    displacement = displacement*np.array(context['spacing'])
    distance = np.linalg.norm(displacement, axis=-1)
    anchors=pairs[...,0,:]
    anchor_B=fields[anchors[...,2],anchors[...,1],anchors[...,0],:3].astype(float)
    denominator=distance*np.linalg.norm(anchor_B,axis=-1)
    cosangle=np.divide(abs(np.sum(displacement*anchor_B,axis=-1)), denominator,
                       out=np.zeros_like(distance),where=denominator>0)
    angle=np.rad2deg(np.arccos(np.clip(cosangle,0,1)))
    rows = []
    for family, indices, ref in groups:
        # The reference is always included, even if all differences happen to vanish.
        ix = list(indices)
        ref_in_group = ix.index(ref)
        labels = [result['labels'][i] for i in ix]
        summary = summarize(delta, tendency[ix], labels, meta['lengths'], meta['directions'], ref_in_group)
        for s, scale in enumerate(meta['scales']):
            for d, direction in enumerate(scale['directions']):
                start = (s*len(meta['directions'])+d)*len(ix)
                block = summary['rows'][start:start+len(ix)]
                for method, row in zip(ix, block):
                    row.update(snapshot=context['snapshot'], time=context['time'], cycle=context['cycle'],
                               nu_coll=context['nu_coll'], scale=scale['scale'], factor=scale['factor'],
                               actual_separation_median=float(np.median(distance[s,d])),
                               actual_angle_to_local_B_median_deg=float(np.median(angle[s,d])),
                               family=family, reference=result['labels'][ref],
                               reference_transfer=block[ref_in_group]['half_structure_function_rate'],
                               paired_transfer_difference=row['transfer_difference_from_reference'])
                    flags = []
                    if direction['below_four_grid_spacings']: flags.append('below 4 cells')
                    if direction['longer_than_half_shortest_box']: flags.append('exceeds half box')
                    if np.any(distance[s,d]==0): flags.append('zero-length pairs')
                    row['resolution_flag'] = ', '.join(flags) or 'resolved'
                    # Compare increments of the signed scalar heat fluxes at the same endpoints.
                    for ch, label in enumerate(('parallel','perp')):
                        cur = np.diff(qs[method,s,d,:,:,ch], axis=-1).ravel()
                        target = np.diff(qs[ref,s,d,:,:,ch], axis=-1).ravel()
                        denom = np.sqrt(np.mean(target**2))
                        row['q_'+label+'_increment_relative_rms'] = float(np.sqrt(np.mean((cur-target)**2))/denom) if denom else None
                    rows.append(row)
    return rows


def trace_case(args, fields, header, dest):
    pars, cfg = header['parameters'], header['config']
    beta = float(pars.get('problem/beta0','nan'))
    nu = float(cfg.get('nu_coll',0))
    dx = np.array(header['spacing'])
    L = min(header['lengths'])
    ds = args.step_cells*min(dx)
    if args.reaches:
        radii = np.rint(np.asarray(args.reaches)*L/ds)*ds
    else:
        # Include physical short-filter candidates as well as box-length tail checks.
        targets = np.r_[args.la,args.ell_nu,L*np.array([.25,.5,1.,2.])]
        radii = np.unique(np.maximum(1,np.rint(targets/ds)))*ds
    if radii[0] <= 0 or np.any(np.diff(radii)<=0):
        raise ValueError('Reaches must be positive and distinct after rounding to the ray step')
    k_band = (np.array(args.fit_band)/L if args.fit_band else
              np.array([0.1/L, np.pi/min(dx)]))
    scales = dict(L_A=args.la, ell_nu=args.ell_nu)
    pairs, meta = build_pairs(fields,dx,scales,args.anchors,args.seed)
    _,_,ids,_ = stencils(pairs, fields.shape[:3])
    ctx = dict(snapshot=header['source'][0],time=header['time'],cycle=header['cycle'],
               nu_coll=nu,spacing=dx.tolist())
    print(f"{len(ids)} unique sampled faces, {round(radii[-1]/ds)} coarse ray steps per arm",flush=True)
    start = perf_counter()
    r = evaluate_faces(fields,ids,dx,cfg,radii,ds,k_band)
    # Same faces, actual physical reaches and fitting coefficients, changing only the ray step.
    fine = evaluate_faces(fields,ids,dx,cfg,radii,ds/2,k_band)
    if (r['labels'] != fine['labels'] or
            not np.array_equal(r['uniform_inverse_scattering_length'],fine['uniform_inverse_scattering_length'])):
        raise RuntimeError('Refinement must hold reaches and uniform scattering fixed')
    n=len(radii)
    samples={}
    groups = [('receiver',range(n+1),n)]
    for j in range(n):
        groups.append(('uniform kernel',[1+n+j,1+2*n+j],1+n+j))
    rows=face_rows(fields,pairs,meta,fine,groups,ctx,samples,'receiver')
    def paired_run(family, indices, q, labels, p=pairs, prefix='step', compare_meta=meta, base=fine):
        merged=dict(q=np.stack(q), labels=labels,bn=base['bn'],bf=base['bf'])
        return face_rows(fields,p,compare_meta,merged,[(family,indices,indices[-1])],ctx,samples,prefix)
    # Step sensitivity of the long reference is part of the same scale report.
    for name, ref in [('receiver',n),('uniform',2*n),('rational',3*n)]:
        rows.extend(paired_run('ray step '+name,[0,1],[r['q'][ref],fine['q'][ref]],
                               [r['labels'][ref]+' (coarse step)',fine['labels'][ref]],prefix='step_'+name))
    optical_note = 'Excluded for zero background scattering' if nu==0 else 'Disabled by --optical-anchors 0'
    if nu>0 and args.optical_anchors:
        optpairs=pairs[:,:,:min(args.optical_anchors,pairs.shape[2])]
        optids=stencils(optpairs, fields.shape[:3])[2]
        optr=args.optical_reach*L
        h=ds/2
        short=optical_faces(fields,optids,dx,cfg,optr/2,h/2)
        long=optical_faces(fields,optids,dx,cfg,optr,h)
        converged=optical_faces(fields,optids,dx,cfg,optr,h/2)
        optical_note=dict(reach=optr, shorter_reach=optr/2, finest_step=converged['actual_ds'],
                          unique_faces=len(optids), path_limiter_fraction=converged['path_limiter_fraction'],
                          path_backup_fraction=converged['path_backup_fraction'], assumptions=converged['assumptions'])
        opt=dict(q=np.stack((converged['q_frozen'],short['q'],long['q'],converged['q'])),
                 bn=converged['bn'],bf=converged['bf'],labels=['receiver, same optical tracer',
                 'optical, half reach','optical, coarse step','optical, full reach'])
        rows+=face_rows(fields,optpairs,meta,opt,[('path scattering',[3,0],0),
                       ('optical reach and step',[1,2,3],3)],ctx,samples,'optical')
    with (dest/'scales.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    np.savez_compressed(dest/'samples.npz',**samples)
    manifest=dict(snapshot=header['source'], shape=header['shape'],beta0=beta,**{k:v for k,v in ctx.items() if k!='snapshot'},
                  parameters=pars, case=args.case, scales_code_units=scales,
                  scale_provenance={'L_A':args.la_source,'ell_nu':args.visc_source}, pairs=meta,
                  ray_reaches_code_units=fine['radii'].tolist(), coarse_ds=ds, fine_ds=ds/2,
                  fit_band_angular_code_units=k_band.tolist(), fit_poles=fine['poles'].tolist(),
                  fit_sampled_relative_error=fine['fit_relative_error'],
                  uniform_inverse_scattering_length=fine['uniform_inverse_scattering_length'].tolist(),
                  optical=optical_note, assumptions=fine['assumptions'],
                  elapsed_analysis_seconds=perf_counter()-start)
    (dest/'metadata.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot',required=True,help='Shared .bin path or quoted glob for all rank files of ONE output')
    p.add_argument('--case',choices=('nu0','collisional'),required=True)
    p.add_argument('--la',type=float,required=True,help='Alfven scale, physical code length')
    p.add_argument('--ell-nu',type=float,required=True,help='Parallel viscous scale, physical code length')
    p.add_argument('--la-source',required=True,help='Source and scale convention; identify matched MHD control')
    p.add_argument('--visc-source',required=True,help='Measurement or explicitly labeled comparison target')
    p.add_argument('--out',type=Path,required=True,help='Fresh output directory')
    p.add_argument('--cache-dir',type=Path,help='Existing scratch directory for temporary field memmap')
    p.add_argument('--anchors',type=int,default=64)
    p.add_argument('--optical-anchors',type=int,default=8,help='Positive-background subset, 0 to disable')
    p.add_argument('--optical-reach',type=float,default=.5,help='Optical reach as a fraction of shortest box side')
    p.add_argument('--reaches',type=float,nargs='+',help='One-sided reaches / shortest box side. Default includes L_A, ell_nu, and .25,.5,1,2 box lengths')
    p.add_argument('--step-cells',type=float,default=.5,help='Coarser of two RK2 ray steps in minimum cell spacings')
    p.add_argument('--fit-band',type=float,nargs=2,metavar=('KMIN_L','KMAX_L'),
                   help='Angular k*L band; default 0.1 to the axial grid Nyquist, L=shortest side')
    p.add_argument('--seed',type=int,default=1729)
    args=p.parse_args()
    if (min(args.la,args.ell_nu,args.step_cells,args.optical_reach)<=0 or
        not np.isfinite([args.la,args.ell_nu,args.step_cells,args.optical_reach,*(args.reaches or []),*(args.fit_band or [])]).all() or
        (args.fit_band and not 0<args.fit_band[0]<args.fit_band[1]) or
        (args.reaches and (len(args.reaches)<2 or args.reaches[0]<=0 or np.any(np.diff(args.reaches)<=0))) or
        args.anchors<2 or args.optical_anchors==1 or args.optical_anchors<0):
        p.error('Need positive finite lengths/steps, >=2 increasing reaches, >=2 anchors, and zero or >=2 optical anchors')
    src=NativeSnapshot(args.snapshot)
    pars=src.metadata['parameters']
    if not np.isclose(float(pars.get('problem/beta0','nan')),16) or pars.get('mhd/cgl_heat_flux')!='landau_fluid':
        p.error('Expected problem/beta0=16 and active landau_fluid; verify the source metadata')
    rate=float(pars.get('mhd/nu_coll','0'))
    if not np.isfinite(rate) or rate<0 or (rate==0)!=(args.case=='nu0'):
        p.error('Wrong collision case for this snapshot')
    field_bytes=int(np.prod(src.shape))*24
    print(f'Staging {field_bytes/2**30:.2f} GiB for the six fields',flush=True)
    if args.cache_dir and shutil.disk_usage(args.cache_dir).free < field_bytes+4096:
        p.error('Insufficient free space for scratch field array')
    if args.out.exists():
        p.error('Output must be a new directory')
    args.out.parent.mkdir(exist_ok=True,parents=True)
    with tempfile.TemporaryDirectory(dir=args.out.parent,prefix='.lf-results-') as tmp:
        dest=Path(tmp)/'result';dest.mkdir()
        if args.cache_dir:
            with tempfile.TemporaryDirectory(dir=args.cache_dir,prefix='lf-fields-') as cache:
                fields,meta=load_fields(src,cache)
                trace_case(args,fields,meta,dest)
                del fields
        else:
            fields,meta=load_fields(src)
            trace_case(args,fields,meta,dest)
        dest.rename(args.out)
    print(args.out/'scales.csv')

if __name__=='__main__':
    main()
