"""Check the IC submodule against the actual GOTHAM C++ potential and cooling tables."""

import io
import logging
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena  # required by the standard regression runner

logger = logging.getLogger('athena' + __name__[7:])
_ROOT = Path(__file__).resolve().parents[3]
_HOST = Path(__file__).resolve().parent / 'cgm_cooling_host'
_EXE = _ROOT / 'tst/build/gotham_ics'


def run():
    _EXE.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++17', '-O2',
                    '-I', str(_HOST), '-I', str(_ROOT / 'src/utils'),
                    '-I', str(_ROOT / 'src/srcterms'),
                    str(_HOST / 'ics.cpp'), '-o', str(_EXE)], check=True)


def analyze():
    sys.path.insert(0, str(_ROOT / 'external/CGM_ICs'))
    import generate_gotham_ics as ic

    config = ic.read_input(_ROOT / 'inputs/gotham/athinput.gotham')
    units = ic.CodeUnits(config)
    potential = ic.GalaxyPotential(config, units)
    Z = config.getfloat('problem', 'metallicity') * (1 - config.getfloat('dust', 'D_floor', fallback=1e-20))
    cooling = ic.Cool.Wiersma_Cooling(Z, 0)
    gravity = ic.G_CGS * units.mass * units.time**2 / units.length**3
    parameters = [gravity] + [config.getfloat('potential', key) for key in
                             ('r_scale', 'rho_scale', 'mass_gal', 'scale_gal',
                              'z_gal', 'r_200', 'rho_mean')]
    rng = np.random.default_rng(42)
    coordinates = rng.normal(size=(2048, 3))
    coordinates *= np.geomspace(1e-4, 400, len(coordinates))[:, None]
    coordinates[0] = 0
    temperature = 10**rng.uniform(3, 8, len(coordinates))
    density = 10**rng.uniform(-7, 3, len(coordinates))
    query = np.column_stack((coordinates, temperature, density, np.full_like(density, Z)))
    stream = io.StringIO()
    stream.write(' '.join(map(str, parameters)) + '\n')
    np.savetxt(stream, query)
    result = subprocess.run([str(_EXE)], input=stream.getvalue(),
                            capture_output=True, text=True, check=True)
    actual = np.loadtxt(io.StringIO(result.stdout))
    expected_phi = potential.phi_cylindrical_kms2(
        np.hypot(coordinates[:, 0], coordinates[:, 1]) * units.kpc,
        coordinates[:, 2] * units.kpc) / units.kms**2
    expected_rate = cooling.lambda_cgs(temperature, density)
    gravity_error = np.max(abs(actual[:, 0] - expected_phi) / np.maximum(abs(expected_phi), 1))
    # The historical C++ table stores rounded log axes; near net-rate zeros use an
    # absolute tolerance instead of dividing by a nearly cancelled rate.
    rate_error = np.max(abs(actual[:, 1] - expected_rate) / (abs(expected_rate) + 1e-23))
    logger.info('IC/C++ gravity max error %.3g; cooling max scaled error %.3g',
                gravity_error, rate_error)
    return gravity_error < 1e-11 and rate_error < 1e-4


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
    raise SystemExit(0 if analyze() else 1)
