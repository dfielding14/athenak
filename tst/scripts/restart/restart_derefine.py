# Regression test for derefinement on restart: restarting an AMR run with a lower
# <mesh_refinement>/num_levels derefines the mesh with total mass unchanged, while a
# restart with the same num_levels leaves the mesh alone.

import logging
from pathlib import Path
import shutil
import subprocess

logger = logging.getLogger('athena' + __name__[7:])

_REPO_DIR = Path(__file__).resolve().parents[3]
_EXE_DIR = _REPO_DIR / 'tst' / 'build' / 'src'
_INPUT = _REPO_DIR / 'inputs' / 'tests' / 'restart_derefine.athinput'
_HST = 'RestartDerefine.hydro.hst'
_OUT_DIR = _EXE_DIR / 'restart_derefine'
_MESSAGE = 'Restart derefined from physical level'


def _athena(args):
    cmd = [str(_EXE_DIR / 'athena')] + args
    logger.debug('Executing: %s', ' '.join(cmd))
    result = subprocess.run(cmd, cwd=str(_EXE_DIR), text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def _masses(path):
    lines = path.read_text().splitlines()
    header = next(line for line in lines if '=mass' in line)
    col = next(int(t.split(']')[0].lstrip('[')) - 1 for t in header.split() if '=mass' in t)
    return [float(line.split()[col]) for line in lines if not line.startswith('#')]


def run(**kwargs):
    for path in (_EXE_DIR / 'rst', _OUT_DIR, _EXE_DIR / 'restart_same', _EXE_DIR / _HST):
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()

    _athena(['-i', str(_INPUT)])
    restart = str(_EXE_DIR / 'rst' / 'RestartDerefine.00000.rst')

    out = _athena(['-r', restart, 'time/nlim=2', '-d', str(_EXE_DIR / 'restart_same')])
    assert _MESSAGE not in out, out

    out = _athena(['-r', restart, 'mesh_refinement/num_levels=1', 'time/nlim=2',
                   '-d', str(_OUT_DIR)])
    assert _MESSAGE + ' 2 to 0' in out, out


def analyze():
    mass0 = _masses(_EXE_DIR / _HST)[0]
    masses = _masses(_OUT_DIR / _HST)
    ok = len(masses) > 0 and all(abs(m/mass0 - 1.0) < 1.0e-12 for m in masses)
    if not ok:
        logger.warning('derefined restart changed total mass: %s vs %s', masses, mass0)
    return ok
