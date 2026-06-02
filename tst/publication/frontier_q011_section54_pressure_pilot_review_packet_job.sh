#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=pic-q011-pressure-review
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

export PYTHONDONTWRITEBYTECODE=1
cd /ccs/home/dfielding/athenak-pic/tst/publication
python3 -B render_q011_section54_pressure_pilot_review_packet.py
