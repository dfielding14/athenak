# Scratch physical-log-ratio prototype (2026-10-01)

This is **design evidence only**, based on commit `7b3345fd2`. It is not a production patch or merge acceptance. The repository worktree was not edited. `q-prototype.patch` contains every scratch change, including the added affine-compression fixture. `source/kokkos` is a symlink to the original checkout's existing Kokkos source.

## Rule implemented

The IAN slot is reinterpreted **globally within this scratch binary** as $Q=\rho q$, where $q=\ln(p_\perp/p_\parallel)$. Primitive/conserved conversion, pressure recovery, HLLE and LLF anisotropy fluxes, and wall re-encoding use this representation. The flux is $F_Q=F_\rho q_{\rm upwind}$. The required equation is

$$
\partial_tQ+\nabla\cdot(Q\boldsymbol v)
=\rho(3\boldsymbol{bb}:\nabla\boldsymbol v-\nabla\cdot\boldsymbol v).
$$

The prototype uses centered cell velocity gradients in 1D, with the same source in the real RK update and FOFC's predicted update. Below `bfloor`, C2P/P2C and face-state conversion set $q=0$ at fixed internal energy; the **entire source** is zero there. Above the floor, pressure recovery uses $q$ directly, and the CGL source resumes. No magnetic reference appears in $q$ encoding or decoding.

Mandatory fluid-firehose wall projection remains enabled. Total-energy and magnetic flux/update algorithms were not altered.

## Contact experiment

`run_prototype.py` runs dc/plm reconstruction, velocities $+10,-10,+1,0$, and `bfloor` values $10^{-14},10^{-10},10^{-6}$: 24 cases, 50 cycles each. Input weak field is always $10^{-12}$, so the $10^{-14}$ case is magnetized everywhere. Each cycle 0 through 50 is checked (Athena writes a duplicate final snapshot).

- All 24 cases pass the unchanged pressure-ratio interval [0.5, 2] with positive pressures.
- Across all snapshots, the actual ratio range is approximately [0.99410, 1.02520].
- Supersonic PLM range is approximately [0.99630, 1.00515].
- There are no density/energy/temperature floor or FOFC events.
- Boundary primitive states remain exactly unchanged.
- The contact has constant initial boundary energy-flux difference $0.25|v|$ and domain length 1. The residual $E_{\rm final}-(E_{\rm initial}+0.25|v|t)$ is at most $7.11\times10^{-15}$ in absolute domain-integrated energy.

Full per-cycle extrema, commands, logs, primitive/conserved-field tables, and histories are retained in `runs/` and summarized in `results.json` / `results.txt`.

## Smooth source discriminator

`run_affine.py` uses an initially uniform plasma with $v_x=a_0x$, $a_0=-0.2$, on [-50,50], evolved to $t=0.2$. The central region |x|<10 is isolated from boundary effects. In the homogeneous exact solution, $\rho=1/(1+a_0t)$, so magnetized CGL predicts

$$
q_\parallel=-2\ln\rho=2\ln(0.96),\qquad
q_\perp=\ln\rho=-\ln(0.96).
$$

Strong uniform fields of magnitude 10 keep these small physical anisotropies away from the mandatory wall. For nx = 128, 256, 512, maximum q errors are:

| Direction | 128 | 256 | 512 |
| --- | ---: | ---: | ---: |
| Parallel | 1.119e-6 | 2.827e-7 | 7.157e-8 |
| Perpendicular | 5.563e-7 | 1.400e-7 | 3.528e-8 |

Both show approximately second-order convergence. Repeating compression with a subfloor transverse field keeps $q=0$ exactly at every resolution. All pressures remain positive, and no floor or FOFC events occur. Results are retained in `affine/`, `affine-results.json`, and `affine-results.txt`.

## Reproduction

`apply_prototype.py` and `add_affine_fixture.py` reconstruct the edits in a fresh git archive of the recorded head. Alternatively apply `q-prototype.patch` to that archive. Build with:

```
cmake -S source -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
python3 run_prototype.py
python3 run_affine.py
```

The run scripts require NumPy, use their containing directory as the experiment root, and contain assertions. Original runs used `/Users/dbf75/.uv/envs/interactive/.venv/bin/python3`. Compiler logs and SHA-256 hashes of scripts, input, patch, results, and final executable are retained.

## Limits

- This prototype uses global Q semantics while retaining old A-oriented internal helper names/comments. Its restarts are not compatible with production A restarts.
- Multidimensional evolution and LF split evolution are rejected at RK update. AMR, MPI, GPU, LF, forcing, viscosity/resistivity, passive CGL, shock heating partitions, and production restart migration are not validated.
- FOFC's predictor and LLF flux were updated consistently, but these experiments never trigger actual FOFC fallback. Forced-fallback tests remain required.
- Centered source gradients suffice for this smooth-source experiment; their behavior at shocks and interaction with reconstruction, AMR, and limiter order require separate design/validation.
- This establishes a viable physical encoding and transition rule. It does not establish that the production architecture should change globally to Q, or that a temporary hyperbolic representation is preferable.
