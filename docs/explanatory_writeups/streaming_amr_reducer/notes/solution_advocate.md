# Solution Advocate

## Strongest Fair Case

The requested products are weighted histograms, so the direct mathematical
operation is a sum over native AMR leaf cells:

`H_p(b) = sum_leaves dV_c w_p(q_c) 1[B_p(x_c,q_c) = b]`.

That sum is naturally partitioned across shards and MPI ranks. It needs only a
bounded input chunk and fixed-size histogram storage. A global dense cube is
not part of the desired quantity and would add interpolation choices plus an
impractical intermediate. At 8 pc across a 400 kpc domain, even eight float32
fields on a uniform grid would be about 4 PB decimal.

The implementation evidence supports the architecture:

- the independent oracle and release-gate suite pass;
- CPU/GPU and one-/two-rank cases agree;
- a complete real shard and an eight-shard/four-GPU run complete;
- same-weight products close to relative spreads of `0` to `6.42e-13`;
- archived controls and the preserved `output29` marginal agree closely;
- memory is bounded by `--chunk-blocks` and histogram count, not snapshot size.

## Claim Boundary

This is strong evidence for the architecture and original repair definitions.
It is not yet evidence that the current global-atomic kernel is fast enough on
MI250X, that all science-product conventions are externally validated, or that
the current logical-key and volume checks prove arbitrary physical coverage.
