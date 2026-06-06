#!/bin/bash
# Directive-only template. The installed trusted trampoline executes the launch contract.
# Q043 uses normal QoS because the exact registered matrix is production-length work.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=pic-q043-registered-execution
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log
