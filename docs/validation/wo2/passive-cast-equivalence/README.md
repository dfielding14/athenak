# Passive invariant cast equivalence

The two complete-expression `static_cast<Real>` wrappers in
`src/eos/cgl_passive.hpp:19-20` produce identical double-precision CPU machine
code and HIP device instructions and constants under the final build flags.
They fix aggregate initialization narrowing when the existing focused helper
checks use `Real=float`. This artifact proves double-precision code equivalence;
it does not claim support for a full float application.

The source comparison asserts that these two wrappers are the only header
changes. Original and changed headers are retained under `cpu/` and `hip/`.
`probe.cpp` exports noinline wrappers around `PassiveSpecificInvariants` and
`PassiveEncode`; the HIP probe also calls both from an exported device kernel.
The compile commands use the exact final CMake `flags.make` settings, with only
one extra header-search directory selecting the old/new header. Both variants
use the same source path and header-search path. The CPU and HIP manifests record
all compiler commands, flags hashes, header hashes and successful return codes.
No application was run for this proof.

| Artifact compared | Result | Bytes | SHA256 of both versions |
| --- | --- | ---: | --- |
| Complete CPU object | Raw bytes identical | 5904 | `eab5aefd0873df5bade5a3f65945135e441bd6fd53ba5a3ee40fccd152cf0fb1` |
| HIP device `.text` | Raw bytes identical | 6868 | `01dc37dc30d09821754155aa70e90336cd9ee0c310832002c2c2d1d092716ec5` |
| HIP device `.rodata` | Raw bytes identical | 32912 | `19a526ad549807ff70d4632ba9fc515ae0296349e453a8d6ff602e9b87354c38` |

Complete emitted CPU IR matches after normalizing only source columns on the
two edited header lines. Complete HIP device IR matches after normalizing only
two occurrences of the compiler-generated `__hip_cuid_*` identifier. Function
bodies and arithmetic attributes match. Raw IR diffs are retained; the original
compile manifests' `passed:false` records that raw metadata differed, rather
than a compiler failure. The separate [equivalence audit](equivalence-audit.json)
passes all machine-code and narrowly normalized IR checks.

Evidence: [CPU compile manifest](cpu/manifest.json),
[HIP compile manifest](hip/manifest.json), [CPU raw IR diff](cpu/ir.diff),
[HIP raw IR diff](hip/ir.diff), [audit source](audit_proof.py),
[probe source](probe.cpp). The raw `.o`, `.ll`, `.hsaco`, extracted sections and
disassemblies are retained alongside them.

Reproduction in the original WO2 scratch directory:

```sh
bash run.sh cpu > cpu-build.log 2>&1
bash run.sh hip > hip-build.log 2>&1
bash extract_device_sections.sh
```

The environment is CCE 20.0.0, ROCm 6.4.2, gfx90a, and the exact frozen final
application compile flags. The final CPU `-mno-daz-ftz` option is a link option;
it is correctly absent from the CPU compile flags. HIP uses the final
`-fno-cray -mno-daz-ftz` flags. No source, application binary, or original header
was modified by this proof.
