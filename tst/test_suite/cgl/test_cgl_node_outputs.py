"""A legacy shared restart can switch its output layout without changing the trajectory."""
import os
import re
from pathlib import Path
import subprocess

import numpy as np

from test_suite.cgl import passive_acceptance


def test_node_restart_round_trip(tmp_path):
    exe = Path(os.environ.get('ATHENAK_TEST_EXE', 'athena')).resolve()
    deck = tmp_path/'initial.athinput'
    deck.write_text(passive_acceptance.input_text(passive=False, active_cgl=True,
                                                  forcing=True, nx=16, nlim=1, tlim=1))

    def run(name, *args):
        folder = tmp_path/name
        folder.mkdir()
        subprocess.run([str(exe), '-d', str(folder), *args], check=True, timeout=120)
        return max((folder/'rst').glob('*.rst'))

    initial = run('initial', '-i', str(deck))
    with initial.open('rb') as pin:
        header = pin.read(32768).split(b'<par_end>')[0].decode()
    rst_block = next(name for name, text in re.findall(r'<(output\d+)>([^<]*)', header)
                     if re.search(r'file_type\s*=\s*rst\b', text))
    shared = run('shared', '-r', str(initial), 'time/nlim=2')
    sharded = run('node', '-r', str(initial), 'time/nlim=2', rst_block+'/single_file_per_node=true')
    shards = list(sharded.parent.glob('node_*/'+sharded.name))
    assert len(shards) == 1 and 0 < sharded.stat().st_size < shards[0].stat().st_size
    reference = run('resume_shared', '-r', str(shared), 'time/nlim=3')
    result = run('resume_node', '-r', str(sharded), 'time/nlim=3', rst_block+'/single_file_per_node=false')
    ref, actual = [passive_acceptance.decode_restart(p, False, active_cgl=True) for p in (reference, result)]
    np.testing.assert_array_equal(ref['active'], actual['active'])
    for a,b in zip(ref['faces'],actual['faces']):
        np.testing.assert_array_equal(a,b)
