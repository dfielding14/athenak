import os
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / 'tst'))
sys.path.insert(0, str(root / 'vis/python'))
binary = Path(sys.argv[1]).resolve()
work = Path(sys.argv[2]).resolve()
os.chdir(root / 'tst')
import pytest
from test_suite.cgl import test_cgl_amr_gpu as amr
amr.INPUT_ROOT = str(root / 'inputs/tests')
work.mkdir(parents=True, exist_ok=True)
link = work / 'athena'
if not link.exists():
    link.symlink_to(binary)
os.chdir(work)
tests = [str(root / 'tst/test_suite/cgl/test_cgl_amr_walls_cpu.py')]
if len(sys.argv) > 3:
    tests += [str(root / 'tst/test_suite/cgl/test_cgl_amr_gpu.py') + '::test_cgl_lf_amr_3d_churn_gpu',
              str(root / 'tst/test_suite/cgl/test_cgl_amr_gpu.py') + '::test_cgl_lf_amr_primitive_churn_gpu',
              str(root / 'tst/test_suite/cgl/test_cgl_c2p_pressure_floor_cpu.py')]
raise SystemExit(pytest.main(['-q', *tests]))
