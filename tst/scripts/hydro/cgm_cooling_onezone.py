# Regression test for the exact operator-split CGM cooling integrator
# (cgm_cooling::IntegrateTemperature): builds cgm_cooling_host/onezone.cpp on the host
# with a shim athena.hpp (no Kokkos) and checks one-zone temperatures against
# reference values. Needs only a C++17 compiler ($CXX, default c++).

import logging
import os
from pathlib import Path
import subprocess

logger = logging.getLogger('athena' + __name__[7:])

_HOST_DIR = Path(__file__).resolve().parent / 'cgm_cooling_host'
_SRC_DIR = Path(__file__).resolve().parents[3] / 'src' / 'srcterms'
_EXE = Path(__file__).resolve().parents[2] / 'build' / 'cgm_cooling_onezone'


def run():
    _EXE.parent.mkdir(parents=True, exist_ok=True)
    cmd = [os.environ.get('CXX', 'c++'), '-std=c++17', '-O2',
           '-I', str(_HOST_DIR), '-I', str(_SRC_DIR),
           str(_HOST_DIR / 'onezone.cpp'), '-o', str(_EXE)]
    logger.debug('Executing: ' + ' '.join(cmd))
    subprocess.run(cmd, check=True)


def analyze():
    result = subprocess.run([str(_EXE)], capture_output=True, text=True)
    logger.info(result.stdout)
    if result.returncode != 0:
        logger.warning('one-zone cooling test: %d case(s) off by > 0.5%%',
                       result.returncode)
    return result.returncode == 0
