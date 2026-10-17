"""Pressure cross spectra supplement; existing turbulence caches are untouched.

Use the turbulence manifest's tasks/cases, with new cache destinations.
Compute MANIFEST, plot cached results with MANIFEST --plot, or --self-check.
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

from analyze_cgl_turbulence import (WORKERS, assemble, bin_power, fourier_grid,
                                    metadata, read_binary)
from plot_cgl_turbulence import time_weights


def cross_spectra(pressure, magnetic_pressure, grid):
    p = pressure - np.mean(pressure, dtype=np.float64)
    m = magnetic_pressure - np.mean(magnetic_pressure, dtype=np.float64)
    phat = fft.rfftn(p, norm='forward', workers=WORKERS)
    mhat = fft.rfftn(m, norm='forward', workers=WORKERS)
    cross = (phat.real*mhat.real + phat.imag*mhat.imag)*grid['weights'][None, None, :]
    curves = dict(P=bin_power(phat, grid), M=bin_power(mhat, grid),
                  C=np.bincount(grid['shells'].ravel(), weights=cross.ravel(),
                                minlength=grid['bins'])/grid['dk'],
                  residual=bin_power(phat+mhat, grid))
    np.testing.assert_allclose(curves['residual'], curves['P']+curves['M']+2*curves['C'],
                               rtol=4e-5, atol=1e-7*max(curves['P'].max(), curves['M'].max()))
    covariance = np.mean(p*m, dtype=np.float64)
    variances = [np.mean(q*q, dtype=np.float64) for q in (p, m)]
    np.testing.assert_allclose([curves[key].sum()*grid['dk'] for key in ('P', 'M')],
                               variances, rtol=4e-5, atol=1e-12)
    scale = np.sqrt(variances[0]*variances[1])
    assert abs(curves['C'].sum()*grid['dk']-covariance) <= 4e-5*scale+1e-12
    return curves


def metrics(curves, band):
    p, m, c, residual = (float(curves[key][band].sum()) for key in ('P', 'M', 'C', 'residual'))
    assert p > 0 and m > 0
    correlation = c/np.sqrt(p*m)
    assert abs(correlation) <= 1+4e-5
    return dict(correlation=correlation, pressure_to_magnetic_power=p/m,
                residual_to_uncorrelated_power=residual/(p+m))


def compute(manifest):
    dependencies = [Path(__file__), Path(__file__).with_name('analyze_cgl_turbulence.py'),
                    Path.home()/'athenak-cgl/vis/python/bin_convert.py']
    codehash = hashlib.sha256(b''.join(path.read_bytes() for path in dependencies)).hexdigest()
    rank, ranks = (int(os.environ.get(key, default)) for key, default in
                   (('SLURM_PROCID', '0'), ('SLURM_NTASKS', '1')))
    for task in manifest['tasks'][rank::ranks]:
        source, output = Path(task['source']), Path(task['cache'])
        output.parent.mkdir(parents=True, exist_ok=True)
        info = output.with_suffix('.json')
        before = source.stat()
        provenance = dict(source=str(source), time=float(task['time']), size=before.st_size, mtime_ns=before.st_mtime_ns,
                          analyzer_sha256=codehash)
        if output.exists() and info.exists() and json.loads(info.read_text())['provenance'] == provenance:
            print('Cached', output, flush=True)
            continue
        start = time.monotonic()
        meta = metadata(source)
        assert before.st_size/(meta['nvars']*meta['varbytes'])*4*70 < 350*2**30
        raw = read_binary(str(source))
        fields, lengths = assemble(raw)
        del raw
        pressure = fields['p_perp']
        magnetic = sum(fields[name]**2 for name in ('bcc1', 'bcc2', 'bcc3'))/2
        del fields
        gc.collect()
        grid = fourier_grid(pressure.shape, lengths)
        curves = cross_spectra(pressure, magnetic, grid)
        curves.update(k=grid['k'], dk=grid['dk'], time=task['time'], actual_time=meta['time'])
        stats = dict(provenance=provenance, time=task['time'], actual_time=meta['time'],
                     k4_32=metrics(curves, (grid['k'] >= 4) & (grid['k'] <= 32)),
                     full=metrics(curves, grid['k'] > 0),
                     elapsed_seconds=time.monotonic()-start,
                     process_peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20)
        after = source.stat()
        assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
        temporary = output.with_suffix('.tmp.npz')
        np.savez(temporary, **curves)
        temporary.replace(output)
        temporary = info.with_suffix('.tmp.json')
        temporary.write_text(json.dumps(stats, indent=2, allow_nan=False)+'\n')
        temporary.replace(info)
        print('Completed', output, stats['elapsed_seconds'], 's', flush=True)
        del pressure, magnetic, grid, curves
        gc.collect()


def plot(manifest, directory):
    assert 1 <= len(manifest['cases']) <= 3, 'Expected one to three cases'
    assert len({case['name'] for case in manifest['cases']}) == len(manifest['cases'])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path.home()/'dbfplot/src'))
    import dbfplot
    dbfplot.use(latex=False, cycle='dark2')
    fig, axes = dbfplot.subplots(1, 2, figsize=(10, 4), layout='constrained')
    summary, common_times = {}, None
    for case, color, style in zip(manifest['cases'], dbfplot.colors('dark2', len(manifest['cases'])), ('-', '--', '-.')):
        snapshots = [dict(np.load(filename, allow_pickle=False)) for filename in case['files']]
        snapshots.sort(key=lambda item: float(item['time']))
        times = np.array([float(item['time']) for item in snapshots])
        if common_times is None:
            common_times = times
        np.testing.assert_array_equal(times, common_times)
        k = snapshots[0]['k']
        for item in snapshots:
            np.testing.assert_array_equal(item['k'], k)
        weights = time_weights(times)
        averaged = {key: np.sum(weights[:, None]*[item[key] for item in snapshots], axis=0)
                    for key in ('P', 'M', 'C', 'residual')}
        band = (k >= 4) & (k <= 32)
        result = dict(times=times.tolist(), weights=weights.tolist(),
                      actual_times=[float(item['actual_time']) for item in snapshots],
                      k4_32=metrics(averaged, band), full=metrics(averaged, k > 0))
        temporal = np.array([metrics(item, band)['correlation'] for item in snapshots])
        mean = float(weights @ temporal)
        result['snapshot_band_correlation'] = dict(mean=mean,
            temporal_std=float(np.sqrt(weights @ (temporal-mean)**2)), min=float(temporal.min()),
            max=float(temporal.max()), values=temporal.tolist())
        for label, keep in [('early', times <= 10), ('late', times >= 10)]:
            selected = np.flatnonzero(keep)
            if len(selected) > 1:
                w = time_weights(times[keep])
                mean_curves = {key: np.sum(w[:, None]*[snapshots[i][key] for i in selected], axis=0)
                               for key in averaged}
                result[label] = dict(window=times[keep][[0, -1]].tolist(), **metrics(mean_curves, band))
        summary[case['name']] = result
        visible = (k > 0) & (k <= manifest.get('kmax', 80)) & (averaged['P']*averaged['M'] > 0)
        axes[0].semilogx(k[visible], averaged['C'][visible]/np.sqrt(averaged['P'][visible]*averaged['M'][visible]),
                         color=color, linestyle=style, label=case['label'])
        axes[1].loglog(k[visible], averaged['residual'][visible]/(averaged['P'][visible]+averaged['M'][visible]),
                       color=color, linestyle=style, label=case['label'])
    for ax in axes:
        ax.axvspan(4, 32, color='0.5', alpha=.10, linewidth=0)
        ax.set(xlabel=r'$kL_x/(2\pi)$', xlim=(float(k[k > 0].min()), manifest.get('kmax', 80)))
        dbfplot.suppress_minor_labels(ax)
    axes[0].set(ylabel=r'$C_{p_\perp,B^2/2}(k)/\sqrt{E_{p_\perp}(k)E_{B^2/2}(k)}$', ylim=(-1.05, 1.05))
    axes[0].axhline(0, color='0.6', linestyle=':', linewidth=1)
    axes[0].legend(frameon=False)
    axes[1].set(ylabel=r'$E_{p_\perp+B^2/2}(k)/[E_{p_\perp}(k)+E_{B^2/2}(k)]$')
    axes[1].axhline(1, color='0.6', linestyle=':', linewidth=1)
    axes[0].text(.97, .97, f'Time average: {common_times[0]:g}–{common_times[-1]:g}',
                 transform=axes[0].transAxes, ha='right', va='top', fontsize=9)
    dbfplot.savefig(fig, directory/'07_pressure_balance', profile='paper', strict=True, ignore=('figure.size',))
    plt.close(fig)
    (directory/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')


def self_check():
    shape, lengths = (64, 32, 32), (1., 1., 2.)
    grid = fourier_grid(shape, lengths)
    z = np.arange(64)[:, None, None]/32
    wave = np.broadcast_to(np.sin(9*np.pi*z), shape).astype(np.float32)
    quadrature = np.broadcast_to(np.cos(9*np.pi*z), shape).astype(np.float32)
    band = (grid['k'] >= 4) & (grid['k'] <= 32)
    for magnetic, correlation, residual in ((-wave, -1, 0), (wave, 1, 2), (quadrature, 0, 1)):
        curves = cross_spectra(3+wave, 4+magnetic, grid)
        values = metrics(curves, band)
        assert np.isclose(values['correlation'], correlation, atol=2e-6)
        assert np.isclose(values['residual_to_uncorrelated_power'], residual, atol=2e-6)
        assert np.isclose(curves['P'][9]*grid['dk'], .5, rtol=3e-6)
    np.testing.assert_allclose(time_weights([6, 7, 10]), [1/8, 1/2, 3/8])
    print('Anti-phase, in-phase, quadrature, rectangular shell, Parseval and time weights passed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', nargs='?', type=Path)
    parser.add_argument('--self-check', action='store_true')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.self_check:
        self_check()
    else:
        manifest = json.loads(args.manifest.read_text())
        if args.plot:
            plot(manifest, args.manifest.parent)
        else:
            compute(manifest)
