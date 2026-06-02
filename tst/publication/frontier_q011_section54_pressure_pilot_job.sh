#!/bin/bash
# Directive-only template. The installed trusted trampoline executes the launch contract.
# Shared by the four separately authorized Q-011 Section 5.4 pressure pilots.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:15:00
#SBATCH --job-name=pic-q011-section54-pressure-pilot
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log
