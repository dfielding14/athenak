# Solution Advocate

## Strongest Fair Case

The remaining uncertainties are common-mode. A wrong HIP build, incomplete
shard set, global AMR overlap, bad 1024-rank reduction, or unusable MI250X
throughput would not create 161 independent failures. It would repeat one
failure 161 times.

The completed tests establish the equations, formats, partial real-data path,
GPU execution, and small MPI reduction. They cannot establish complete shard
identity, global leaf uniqueness, ancestor exclusion, domain-volume closure,
or production-hardware throughput. One exact controlled snapshot is therefore
the cheapest test that reaches the remaining failure modes.

The chosen golden snapshot also has archived controls. A successful run must
match intact `output31`, match the preserved radial marginal of `output29`,
disagree with its broken angular distribution, activate global closure checks,
and expose the relative costs of `original`, `science`, and `all` products.

## Claim Boundary

The strategy is strongly supported, but it has not yet been demonstrated on
Frontier or enforced by the current submitter. One golden success also does
not remove per-snapshot validation, the need for a first 256-node canary, or
the need to freeze the exact executable that passed.
