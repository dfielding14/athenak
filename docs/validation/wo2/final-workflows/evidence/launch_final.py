#!/usr/bin/env python3
"""Run only the frozen final Athena binary in a coordinated GPU Slurm step."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time

binary = Path(os.environ['WO2_FINAL_BINARY']).resolve()
expected = os.environ['WO2_FINAL_SHA256']
if hashlib.sha256(binary.read_bytes()).hexdigest() != expected:
    raise SystemExit('Final binary changed after preflight')
args = sys.argv[1:]
if Path(sys.argv[0]).name == 'mpirun':
    if len(args) < 3 or args[0] not in ('-np', '-n'):
        raise SystemExit('Final mpirun adapter expects -np/-n RANKS COMMAND')
    ranks = int(args[1])
    if Path(args[2]).resolve() != Path(__file__).resolve():
        raise SystemExit('Final mpirun adapter accepts only the retained launcher')
    command = args[2:]
else:
    if os.environ.get('SLURM_STEP_ID', '') not in ('', 'batch', 'extern'):
        os.execv(str(binary), [str(binary), *args])
    ranks = 1
    command = [str(binary), *args]
if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Final GPU acceptance requires a coordinated Slurm allocation')
launch = ['srun', '--exact', '-N', '1', '-n', str(ranks),
          '--ntasks-per-node', str(ranks), '--threads-per-core=1',
          '--cpu-bind=threads', '-c', '7', '--gpus-per-task=1', '--gpu-bind=closest']
log = Path(os.environ['WO2_FINAL_LAUNCH_LOG'])
with log.open('a') as stream:
    stream.write(json.dumps({'utc_epoch':time.time(), 'cwd':str(Path.cwd()),
        'binary':str(binary), 'sha256':expected, 'command':[*launch, *command],
        'allocation':os.environ['SLURM_JOB_ID']}) + '\n')
os.execvp('srun', [*launch, *command])
