#!/bin/sh
set -eu
export PYTHONDONTWRITEBYTECODE=1
/Users/dbf75/.uv/envs/interactive/.venv/bin/python3 /Users/dbf75/.codex/worktrees/c43e/athenak-DF/benchmarks/reconnection/run_local.py --binary /Users/dbf75/.codex/worktrees/c43e/athenak-DF/build-reconnection-release/src/athena --output /Users/dbf75/.codex/worktrees/c43e/athenak-DF/build-reconnection-release/campaign/cl-cpd2-t5 --restart /Users/dbf75/.codex/worktrees/c43e/athenak-DF/build-reconnection-release/campaign/cl-cpd2-t1/rst/resistive_harris.00002.rst --model current_limited --cells-per-di 2 --tlim 5 --cfl .4 --sts-ratio 32 --output-dt .1 --history-dt .01 --ranks 3 --wall-limit 00:15:00
