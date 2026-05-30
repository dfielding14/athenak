#!/bin/bash
# Directive-only template. The installed trusted trampoline executes the launch contract.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:15:00
#SBATCH --job-name=pic-f1-structured-gpu-gyro
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log
