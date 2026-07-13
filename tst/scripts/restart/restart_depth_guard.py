import logging
from pathlib import Path
import shutil
import subprocess

import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_REPO_DIR = Path(__file__).resolve().parents[3]
_EXE_DIR = _REPO_DIR / 'tst' / 'build' / 'src'
_INPUT = _REPO_DIR / 'inputs' / 'tests' / 'restart_depth_guard.athinput'


def _run(args, cwd, check=True):
    cmd = [str(_EXE_DIR / 'athena')] + args
    return subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, check=check)


def run(**kwargs):
    logger.debug('Running test %s', __name__)
    rst_dir = _EXE_DIR / 'rst'
    same_depth_dir = _EXE_DIR / 'restart_same_depth'
    for path in (rst_dir, same_depth_dir):
        if path.is_dir():
            shutil.rmtree(path)

    _run(['-i', str(_INPUT)], _EXE_DIR)
    restart = rst_dir / 'RestartDepthGuard.00000.rst'
    assert restart.is_file(), 'missing source restart for depth-guard test'

    _run(['-r', str(restart), 'time/nlim=0', 'time/tlim=0.0',
          '-d', str(same_depth_dir)], _EXE_DIR)

    rejected = _run(['-r', str(restart), 'mesh_refinement/num_levels=1',
                     'time/nlim=0', 'time/tlim=0.0'], _EXE_DIR, check=False)
    assert rejected.returncode != 0, 'ordinary over-deep restart unexpectedly succeeded'
    output = rejected.stdout + rejected.stderr
    assert 'A restart conversion is required' in output
    assert 'logical level 2' in output


def analyze():
    return True
