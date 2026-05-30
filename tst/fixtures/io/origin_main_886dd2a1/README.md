# Origin Main IO Baseline Fixtures

These files are immutable compatibility fixtures generated from the clean baseline:

- source revision: `origin/main` at `886dd2a1437e45a3a30b3eeebf2adfa838328f73`;
- serial build: `cmake -S /tmp/athenak-io-baseline-src -B /tmp/athenak-io-baseline-build -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=OFF`;
- MPI build: `cmake -S /tmp/athenak-io-baseline-src -B /tmp/athenak-io-baseline-build-mpi -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=ON`;
- serial producer: `producer/origin_main_legacy_shared.athinput`;
- serial two-dimensional PDF producer: `producer/origin_main_legacy_pdf_2d.athinput`;
- two-rank producer: `producer/origin_main_legacy_per_rank.athinput`.

The producer commands were:

```sh
/tmp/athenak-io-baseline-build/src/athena -i producer/origin_main_legacy_shared.athinput -d /tmp/athenak-io-baseline-shared
/tmp/athenak-io-baseline-build/src/athena -i producer/origin_main_legacy_pdf_2d.athinput -d /tmp/athenak-io-baseline-pdf2d
mpirun -np 2 /tmp/athenak-io-baseline-build-mpi/src/athena -i producer/origin_main_legacy_per_rank.athinput -d /tmp/athenak-io-baseline-rank
```

The artifacts establish the legacy contract for:

| Directory | Baseline behavior |
| --- | --- |
| `bin/shared/` | Existing serial/shared binary output at initial and terminal output times. |
| `bin/per_rank/` | Existing two-rank `single_file_per_rank` binary partitioning. |
| `cbin/shared/` | Existing serial/shared coarsened binary output. |
| `cbin/per_rank/` | Existing two-rank `single_file_per_rank` coarsened binary partitioning. |
| `pdf/legacy_1d/` | Existing text-form one-dimensional PDF values and companion `.bins.pdf`. |
| `pdf/legacy_2d/` | Existing text-form two-dimensional PDF values and companion `.bins.pdf`. |
| `rst/shared/` | Existing serial/shared restart output. |
| `rst/per_rank/` | Existing two-rank `single_file_per_rank` restart partitioning. |

`MANIFEST.tsv` records relative path, byte size, and SHA-256 digest for every
producer deck and generated output; `SHA256SUMS` permits direct digest checking.
New N-dimensional PDF, spherical-slice, and per-node formats are intentionally
absent: they must be tested against their new reader contracts and cross-layout
numerical comparisons, not misrepresented as baseline formats.
