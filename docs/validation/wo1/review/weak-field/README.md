# T-B4 multi-cycle weak-field transport: unresolved merge blocker

The retained B4 correction prevents the immediate first-cycle contamination but
fails the requested multi-cycle grid check. Both flow directions fail identically.
The regression is marked `xfail(strict=True)` with a merge-blocking explanation;
its pressure-ratio bounds remain [0.5, 2]. An unexpected pass fails CI until the
marker is removed after reviewing a fix.

The 128-cell initial state has constant density 1 and velocity ±10, zero normal
field, and a transverse field jumping from $10^{-12}$ to 1. Both sides start
isotropic. Gas pressures 1.5 and 1 give the same total pressure $p+B^2/2=1.5$.
`bfloor=1e-10`; collisions and kinetic limiters are off. The contact stays away
from the outflow boundaries over 50 RK2 cycles. The exact advected discontinuity
retains its isotropic states. This is an application-level test, not the original
single-cell flux update.

| Flux headers | Reconstruction | Maximum ratio at cycle 1 | Maximum ratio at cycle 50 |
|---|---|---:|---:|
| Current | PLM | 1 | 2.4630400339e10 |
| Current | Donor cell | 1 | 5.4662083152e5 |
| Pre-B4 | PLM | 1.1098016737e8 | 2.2645066155e12 |
| Pre-B4 | Donor cell | 1.1261634707e8 | 2.2793349357e12 |

The current PLM result first exceeds 2 at cycle 5 (3.00522). At its cycle-50
maximum, $B=4.0462\times10^{-6}>B_{\rm floor}$,
$p_\parallel=9.2342\times10^{-11}$ and $p_\perp=2.27442$. The donor-cell result
shows that this is not solely a high-order reconstruction issue. The precise
flux redesign remains unresolved; no new flux correction is included here.

`results.json` retains every cycle, both directions, all eight cases, executable
hashes, and exact commands. `summary.log` and `pre-b4-build.log` retain the
comparison and compile/link evidence. The pre-B4 control replaces only the HLLE
and LLF flux headers with the parent of task commit `7d1f79557`; it uses the
current application and new problem generator otherwise.

From the repository root, with a current Release CPU Makefiles build in
`tst/build` and NumPy/pytest available:

```sh
python docs/validation/wo1/review/weak-field/build_pre_b4.py . /tmp/wo1-weak-field
python docs/validation/wo1/review/weak-field/b4_compare.py . /tmp/wo1-weak-field
(cd tst/build/src && python -m pytest ../../test_suite/cgl/test_cgl_weak_field_cpu.py --runxfail)
```

The last command deliberately exposes both failing assertions. Without
`--runxfail`, pytest reports two expected failures, not two passing checks.
