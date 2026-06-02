#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=pic-q011-pressure-publish
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PYTHONDONTWRITEBYTECODE=1

cd /ccs/home/dfielding/athenak-pic/tst/publication

python3 -B publish_q011_section54_pressure_pilot_bundle.py \
  "$PIC_ROOT/publication/q011_section54_pressure_pilot_bundle" \
  --receipt-path \
  "$PIC_ROOT/publication/q011_section54_pressure_pilot_bundle_receipt.json" \
  --analysis-result-path \
  "$PIC_ROOT/publication/q011_section54_pressure_pilot_analysis.json" \
  --case-artifact-dir \
  "ps_p0_1p00=$PIC_ROOT/runs/q011_section54_pressure_ps_p0_1p00/a9db403f-a174-45ae-956d-f96185a94db1" \
  --case-descriptor-sha256 \
  "ps_p0_1p00=da322586da4831c7747e3a2f5774adf6b999835dc8636b9d971411e1dc9471b2" \
  --case-artifact-dir \
  "ps_p0_0p05=$PIC_ROOT/runs/q011_section54_pressure_ps_p0_0p05/3eb84dfa-0be7-4a62-b77b-64131c08dc6f" \
  --case-descriptor-sha256 \
  "ps_p0_0p05=978cfc63efc9503b9dd0a043793d759e930f390340b843a39d621c2b859ff85c" \
  --case-artifact-dir \
  "ps_p0_0p10=$PIC_ROOT/runs/q011_section54_pressure_ps_p0_0p10/07bd4d96-ffc8-4446-98ef-562d4d83da0c" \
  --case-descriptor-sha256 \
  "ps_p0_0p10=dba35e7b157c1646444792cfda68963dbfd2062b4887961bfb380e51a113f5d4" \
  --case-artifact-dir \
  "ps_p0_0p20=$PIC_ROOT/runs/q011_section54_pressure_ps_p0_0p20/8b4065a3-b02a-4312-8238-20d72d958bbb" \
  --case-descriptor-sha256 \
  "ps_p0_0p20=ef842e1440cbe1fd51e3ee79af5c838704e782708e2a827ec0bfcc89ca1ae312"

python3 -B publish_q011_section54_pressure_pilot_bundle.py \
  --verify-published-receipt \
  "$PIC_ROOT/publication/q011_section54_pressure_pilot_bundle_receipt.json"
