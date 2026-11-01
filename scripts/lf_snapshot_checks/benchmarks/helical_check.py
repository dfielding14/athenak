"""Offline linear SHD kernel check. Run with python3 helical_check.py."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import RegularGridInterpolator
from scipy.special import sici

root = Path(__file__).resolve().parent
pi = np.pi
amplitude = 1.0
stretch = np.sqrt(1 + amplitude**2)
kparallel = 2*pi/stretch
starts = np.random.default_rng(721).uniform(0, 1, (48, 3))


def temperature(points, mode):
    x, y, z = np.moveaxis(points, -1, 0)
    X = x - amplitude*np.sin(2*pi*z)/(2*pi)
    Y = y + amplitude*np.cos(2*pi*z)/(2*pi)
    return np.cos(2*pi*(X + Y + mode*z))


def trace(points, distance):
    z0 = points[:, 2, None]
    z = z0 + distance[None, :]/stretch
    x = points[:, 0, None] + amplitude*(np.sin(2*pi*z)-np.sin(2*pi*z0))/(2*pi)
    y = points[:, 1, None] - amplitude*(np.cos(2*pi*z)-np.cos(2*pi*z0))/(2*pi)
    return np.stack(np.broadcast_arrays(x, y, z), axis=-1)


def kernel(s, collision):
    if collision == 0:
        return 1/s
    si, ci = sici(collision*s)
    return 1/s-collision*(np.sin(collision*s)*ci-np.cos(collision*s)*(si-pi/2))


def exact_response(reach, collision):
    value = quad(lambda s: kernel(s, collision)*np.sin(kparallel*s), 0,
                 reach, epsabs=1e-11, limit=300)[0]
    return 2*value/pi


records = []
for cells in (32, 64, 128):
    nodes = np.arange(cells)/cells
    grid = np.stack(np.meshgrid(nodes, nodes, nodes, indexing='ij'), axis=-1)
    axes = (np.arange(cells+1)/cells,)*3
    interpolators = {}
    for mode in (0, 1):
        values = temperature(grid, mode)
        values = np.pad(values, ((0, 1),)*3, mode='wrap')
        interpolators[mode] = RegularGridInterpolator(axes, values)
    for reach in (0.125, 0.5, 2.0):
        panels = int(np.ceil(reach*cells))
        knots, weights = np.polynomial.legendre.leggauss(8)
        ds = reach/panels
        distances = ((np.arange(panels)[:, None]+0.5+knots/2)*ds).ravel()
        weights = np.broadcast_to(weights*ds/2, (panels, 8)).ravel()
        forward, backward = trace(starts, distances), trace(starts, -distances)
        for collision in (0.0, 2.0):
            multiplier = exact_response(reach, collision)
            for mode in (0, 1):
                delta = interpolators[mode](forward % 1)-interpolators[mode](backward % 1)
                flux = -np.sum(delta*kernel(distances, collision)*weights, axis=1)/pi
                analytic_delta = temperature(forward, mode)-temperature(backward, mode)
                traced_flux = -np.sum(analytic_delta*kernel(distances, collision)*weights, axis=1)/pi
                # q/(C rho c epsilon)=+response*sin(phase), so transport is down-gradient.
                x, y, z = starts.T
                phase = 2*pi*(x-amplitude*np.sin(2*pi*z)/(2*pi)
                              +y+amplitude*np.cos(2*pi*z)/(2*pi)+mode*z)
                expected = multiplier*np.sin(phase) if mode else np.zeros(len(starts))
                records.append(dict(cells=cells, reach=reach, collision_inverse_length=collision,
                                    mode=mode, finite_response=multiplier,
                                    infinite_response=kparallel/(kparallel+collision),
                                    exact_trace_max_error=float(np.max(abs(traced_flux-expected))),
                                    interpolation_rms_error=float(np.sqrt(np.mean((flux-expected)**2))),
                                    interpolated_flux_rms=float(np.sqrt(np.mean(flux**2))),
                                    expected_flux_rms=float(np.sqrt(np.mean(expected**2)))))

# A direct conservation check of the linear Fourier generator in field coordinates.
z = np.arange(1024)/1024
wave = np.cos(2*pi*z)
rate = kparallel**2/(kparallel+2)
generator = -rate*wave
assert abs(np.mean(generator)) < 1e-13
assert np.mean(wave*generator) < 0
assert max(r['exact_trace_max_error'] for r in records) < 1e-8
for reach in (.125, .5, 2.):
    for collision in (0., 2.):
        chosen = [r for r in records if r['reach']==reach and r['collision_inverse_length']==collision]
        for mode in (0, 1):
            errors = [r['interpolation_rms_error'] for r in chosen if r['mode']==mode]
            assert errors[-1] < errors[0]/8, errors
result = dict(model='Linear uniform-coefficient SHD, helical constant-strength field',
              field_amplitude=amplitude, kparallel=kparallel,
              sampling='48 fixed random starts; periodic trilinear interpolation; 8-point Gauss panels',
              flux_normalization='q/(C rho c epsilon); T=T0+epsilon*cos(phase)',
              conservation_mean=float(np.mean(generator)),
              dissipation_inner_product=float(np.mean(wave*generator)), records=records)
destination = root/'helical_results.json'
destination.write_text(json.dumps(result, indent=2)+'\n')
print(destination)
for row in records:
    if row['cells']==128:
        print(row['reach'], row['collision_inverse_length'], row['mode'],
              round(row['finite_response'], 7), f"{row['interpolation_rms_error']:.3g}")
