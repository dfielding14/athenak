# CR-Shock Storyline Final Report

## Executive Summary
Implemented and validated a laptop-feasible (`np=8`) CR-shock storyline workflow from deck setup through diagnostics, scan orchestration, and figure generation.

Key outcomes:
- Local publication deck and HPC-target deck updated for shock storyline use.
- New deterministic diagnostics extraction pipeline (`shock profile`, `f(x,p)` proxy, `dN/dp`).
- Scan runner upgraded with local/HPC tiers, launcher templating, and JSON/CSV summaries.
- Shock-only 4-panel figure pipeline implemented and reproducibility-checked.
- Final validation reruns completed: smoke, local canonical, and scan subset.

## What Changed (Files + Purpose)
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_CR_SHOCK_STORYLINE_EXECUTION_PLAN.md`
  - Live step-by-step execution tracker, gates, and artifact log.
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_amr_shock_lb_publication_local.athinput`
  - Local/laptop profile tuned for stable long runs on `np=8`; rich output cadence for diagnostics.
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_amr_shock_lb_publication_hpc.athinput`
  - HPC-target high-dynamic-range profile (`384x384x96`, `ppc=16`, `tlim=10.0`, `nlim=4000`).
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/publication/analyze_pic_shock_storyline.py`
  - New diagnostics extractor for shock position series, upstream/downstream windows, `x-p` proxy, and spectra.
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/publication/run_pic_shock_scan.py`
  - Tiered scan matrix runner with `np=8` local default, HPC-scale defaults (`nproc=4096`), launcher templates, and JSON/CSV summaries.
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/publication/plot_pic_shock_storyline.py`
  - New 4-panel storyline figure script (fluid structure, CR phase-space, CR spectrum, AMR/LB context).
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_CR_SHOCK_STORYLINE_FINAL_REPORT.md`
  - This final handoff report.

## Scan Matrix and Outcomes
### Local mini-scan (`np=8`)
Source:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step4/scan_summary.json`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step4/scan_summary.csv`

Executed subset (3 cases; all passed):
- `pic_shock_scan_local_m0p80_ppc2p00_n32` (`runtime_s ~ 1.7`)
- `pic_shock_scan_local_m0p80_ppc2p00_n48` (`runtime_s ~ 4.2`)
- `pic_shock_scan_local_m0p80_ppc2p00_n64` (`runtime_s ~ 8.4`)

### HPC command emission (thousands-of-GPU ready)
Source:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724_hpc_emit/reviews/shock_step4/scan_commands.sh`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724_hpc_emit/reviews/shock_step4/scan_summary.json`

Emission profile:
- `--tier hpc --nproc 4096 --launch-template "srun -n {nproc} --gpus-per-task={gpus_per_rank}" --gpus-per-rank 1`
- Produces scheduler-ready command lines for direct supercomputer execution.

## Canonical Runs and Why
Canonical selection record:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step6/canonical_selection.md`

Selected:
1. Local canonical: `pic_amr_shock_pub_local_step2_np8`
   - Stable on `np=8`, complete diagnostics availability, deterministic pipeline outputs.
2. HPC candidate: `pic_amr_shock_lb_publication_hpc.athinput`
   - Higher resolution/ppc and longer horizon; intended for large-rank GPU campaigns.

## Figure Pack Path
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step5/shock_storyline_multipanel.png`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step5/figure_summary.json`

## Final Validation and Reproducibility
Step 7 validation artifacts:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step7/smoke_rerun_metrics.json` (smoke pass)
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step7/local_canonical_rerun.log` (local canonical pass)
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724_step7_subset/reviews/shock_step4/scan_summary.json` (scan subset pass)
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step7/error_scan.txt` (empty)

Determinism checks:
- Diagnostics checksum diff empty:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step3/checksum_diff.txt`
- Figure checksum diff empty:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step5/figure_sha_diff.txt`

## Residual Risks and Next HPC Actions
Residual risks:
- Adaptive AMR churn can trigger repeated `particle orphaned` messages in longer PIC runs on this laptop-scale workflow.
- Local canonical currently defers AMR churn (`ncycle_check=1000`, `refinement_interval=1000`) for stability.

Next HPC actions:
1. Re-enable active AMR churn in HPC campaign runs and validate no orphan bursts at scale.
2. Use emitted HPC command templates (`nproc=4096`, scheduler launch prefix) for cluster execution.
3. Re-run restart parity and load-balance trend checks under production launcher.
4. Expand scan execution beyond mini-subset to full matrix for publication selection.
