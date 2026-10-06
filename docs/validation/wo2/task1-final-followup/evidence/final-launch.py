#!/usr/bin/env python3
"""Launch retained AthenaK binaries in Slurm steps, including legacy mpirun calls."""
import os
from pathlib import Path
import sys

binary = os.environ['WO2_BINARY']
backend = os.environ.get('WO2_BACKEND', 'hip')
arguments = sys.argv[1:]
if Path(sys.argv[0]).name == 'mpirun':
    if len(arguments) < 3 or arguments[0] not in ('-np', '-n'):
        raise SystemExit('WO2 mpirun adapter expects -np/-n RANKS COMMAND')
    ranks = int(arguments[1])
    command = arguments[2:]
else:
    if os.environ.get('SLURM_STEP_ID', '') not in ('', 'batch', 'extern'):
        os.execv(binary, [binary, *arguments])
    ranks = 1
    command = [binary, *arguments]

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('WO2 application tests require a Slurm allocation')
launch = ['srun', '--exact', '-N', '1', '-n', str(ranks),
          '--ntasks-per-node', str(ranks), '--threads-per-core=1',
          '--cpu-bind=threads', '-c', '7' if backend == 'hip' else '1']
if backend == 'hip':
    launch += ['--gpus-per-task=1', '--gpu-bind=closest']
else:
    launch += ['--gpus-per-node=0', '--overlap']
os.execvp('srun', [*launch, *command])
