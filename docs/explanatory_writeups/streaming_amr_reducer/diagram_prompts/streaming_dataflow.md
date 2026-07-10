# Streaming AMR Data Flow

**Purpose:** Show the complete path from native AMR shards to published PDFs without a dense cube.

**Image-generation prompt:**

> Left-to-right architecture diagram for a high-performance scientific reducer. On the left, many AMR shard files labeled node_00000000 through node_N. In the center, several MPI rank lanes, each assigned different shards. Within each lane show repeated stages: POSIX pread chunk, host buffer, host-to-GPU copy, Kokkos derive fields, Kokkos histogram atomics. On the right, local histograms are combined by MPI reduction, validated, and published as repair/science PDFs plus a manifest. Make clear that full 3D blocks stream through memory and are not assembled into a global dense cube. Include callout: “memory scales with chunk size plus histograms, not full cube.”

**Caption:** Native AMR MeshBlocks stream through independent ranks and GPUs; only histograms are combined globally.

**How to read it:** Follow one bounded chunk. It is read, binned, discarded, and replaced.

**Evidence status:** Schematic, not evidence.
