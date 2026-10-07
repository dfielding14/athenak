# Three-cell mirror tail at t=14.250230767343458

One shared float32 primitive snapshot was read through the existing production analyzer reader in a CPU-only `srun --exact --overlap` step in allocation 5629677. No solver was run and no physics, tolerances, thresholds or snapshot bytes were changed. Snapshot SHA256: `41232a8f3d2000935522d9e1cdfa6f0a0dddc729823ba88ac480428268819e17`.

The audit reproduces exactly 3 cells with `X-envelope > 1` among 1,769,472 cells (volume fraction `1.6954210069444444e-06`). Raw strict `X>1` counts 4,872 cells; most do not exceed the independently propagated float32 rounding envelope. Robust firehose crossings are zero.

| Global (k,j,i) | X−1 | Rounding envelope | Excess beyond envelope | Local scalar beta |
| --- | ---: | ---: | ---: | ---: |
| (93, 10, 41) | 0.000224516823 | 0.000201126811 | 2.3390012e-05 | 2964.84273 |
| (94, 11, 44) | 0.000213621441 | 0.000212454203 | 1.16723775e-06 | 3151.01672 |
| (104, 93, 24) | 4.02660344e-06 | 4.01795962e-06 | 8.64381492e-09 | 62.9136769 |

`cells.json` retains exact rho, both pressures, all three B components, B² and pressure-anisotropy excesses. The first two cells have B²≈0.00474 and 0.00449 and beta≈2965 and 3151; the third has B²≈0.24135 and beta≈62.9. Their absolute `p_perp-p_parallel-B²/2` values are only 5.325e−7, 4.797e−7, 4.859e−7. The rounding bound is particularly wide for the two weak-field cells because X divides the pressure difference by B². These three robust excesses should remain reported; the envelope is not a tolerance that changes the physical threshold.

The retained input has `nu_coll=0`, `limiter_nu_coll=1e10`, `mirror_threshold=1`, and `backup_limiters=false`. In `src/eos/ideal_c2p_mhd.hpp::SingleCollRates_CGLMHD`, positive mirror excess is relaxed as `delta_p = threshold + (delta_p-threshold)/(1+nu_lim*dt)`. A finite rate therefore permits a positive residual; it is not a projection onto X≤1. `SingleCollWalls_CGLMHD` supplies no finite mirror upper wall when backups are disabled. The unconditional physical hard wall is the distinct firehose lower bound X≥−2.

Source order: `MHD::STSPostSweepCGLCollisions` applies full rates once after the post sweep using the cycle duration, then records admissibility. `Driver::Execute` advances time and writes outputs afterward; `sts_merge_half_sweeps=false` in this run. Thus this snapshot samples the post-operator state, not a retained unrelaxed intermediate stage. The measured small positive mirror tails are compatible with finite soft-limiter semantics. This one saved float32 state does not retain pre-relaxation excess or prove that finite relaxation alone accounts for every rounding-scale residual.

This is separate from the old firehose restart bug: there is no violated configured mirror hard wall, no restart decode comparison here, and no suggested threshold relaxation. Clean strict hard-volume/floor diagnostics do not imply zero soft-mirror excess. Likewise `lf_hwproj=0` is not evidence of absent projections; that counter is uninstrumented.

`audit_snapshot.py`, `run.sh`, `run.log`, and `cells.json` retain the calculation, allocated launch, full values and source/snapshot hashes. The initial attempt named the script `inspect.py`, shadowing Python’s standard module; its nonzero-return log is preserved as `run-initial-module-shadow.log`. Renaming the script resolved this harness issue; the successful rerun changed no calculation.
