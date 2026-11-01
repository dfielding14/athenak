"""Simplified constant-speed optical-depth SHD reference on a periodic line.

Positive Krook scattering only. This prescribed-temperature flux check does not
validate a nonlinear CGL closure or evolution across transparent intervals.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline

root = Path(__file__).resolve().parent
pi = np.pi
receiver = .65
base_nu_over_c = 7.0


def flux(values, spacing, inverse_length):
    k = 2*pi*np.fft.fftfreq(len(values), d=spacing)
    multiplier = -1j*k/(abs(k)+inverse_length)
    return np.fft.ifft(multiplier*np.fft.fft(values)).real


records = []
for cells in (4096, 8192, 16384):
    x = np.arange(cells+1)/cells
    # Periodized smooth drive and scattering barrier, separated from receiver.
    temperature = np.sum([np.exp(-((x-.25+j)/.04)**2) for j in (-1, 0, 1)], axis=0)
    bump = np.sum([np.exp(-((x-.45+j)/.035)**2) for j in (-1, 0, 1)], axis=0)
    temperature[-1] = temperature[0]
    for name, coefficient in (('parallel', (3*pi-8)/np.sqrt(8*pi)),
                              ('perpendicular', np.sqrt(2/pi))):
        a0 = coefficient*base_nu_over_c
        frozen = flux(temperature[:-1], 1/cells, a0)
        frozen_at_receiver = CubicSpline(x, np.r_[frozen, frozen[0]], bc_type='periodic')(receiver)
        for height in (0., 1., 9., 99.):
            a = a0*(1+height*bump)
            eta = np.r_[0., cumulative_trapezoid(a, x)]
            optical_period = eta[-1]
            count = 4*cells
            eta_grid = np.arange(count)*optical_period/count
            temperature_eta = CubicSpline(eta, temperature, bc_type='periodic')(eta_grid)
            q_eta = flux(temperature_eta, optical_period/count, 1.)
            q_on_x = CubicSpline(np.r_[eta_grid, optical_period],
                                np.r_[q_eta, q_eta[0]], bc_type='periodic')(eta)
            q_receiver = float(CubicSpline(x, q_on_x, bc_type='periodic')(receiver))
            rows = dict(cells=cells, channel=name, background_inverse_length=float(a0),
                        barrier_peak_factor=1+height, optical_period=float(optical_period),
                        q_receiver=q_receiver, frozen_coefficient_q_receiver=float(frozen_at_receiver),
                        receiver_ratio_to_frozen=float(q_receiver/frozen_at_receiver),
                        receiver_scattering_factor=float(1+height*np.exp(-((receiver-.45)/.035)**2)),
                        uniform_reference_max_error=float(np.max(abs(q_on_x[:-1]-frozen))))
            if height == 0:
                assert rows['uniform_reference_max_error'] < 1e-8
            records.append(rows)

for name in ('parallel', 'perpendicular'):
    values = [r['q_receiver'] for r in records if r['channel']==name and r['barrier_peak_factor']==100]
    assert abs(values[-1]-values[-2])/abs(values[-1]) < 1e-4
    assert np.isfinite(values[-1])
result = dict(model='Constant-speed prescribed-T optical-depth SHD, strictly positive scattering',
              x_period=1., temperature_center=.25, temperature_width=.04,
              scattering_barrier_center=.45, scattering_barrier_width=.035,
              receiver=receiver, flux_normalization='q/(C rho c delta T)',
              limitation='Simplified positive Krook scattering; no nonlinear CGL or transparent intervals',
              records=records)
destination = root/'barrier_results.json'
destination.write_text(json.dumps(result, indent=2)+'\n')
print(destination)
for record in records:
    if record['cells']==16384:
        print(record['channel'], record['barrier_peak_factor'], record['q_receiver'],
              record['receiver_ratio_to_frozen'])
