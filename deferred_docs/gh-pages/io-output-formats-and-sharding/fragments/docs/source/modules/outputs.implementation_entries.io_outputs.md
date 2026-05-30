## Implementation Entry Points

| File | Responsibility |
| --- | --- |
| `src/outputs/outputs.cpp` | Block parsing, writer dispatch, compatibility checks |
| `src/outputs/basetype_output.cpp` | Resolves selected variables and gathers output data |
| `src/outputs/derived_variables.cpp` | Computed variables requested through output names |
| `src/outputs/history.cpp` | Integrated diagnostics |
| `src/outputs/binary.cpp` | Native shared, rank-sharded, and node-sharded mesh binary writing |
| `src/outputs/coarsened_binary.cpp` | Uniform 3D active-zone coarsened-binary writing |
| `src/outputs/restart.cpp` | Checkpoint writing and read-compatible state |
| `src/outputs/pdf.cpp` | One- through four-dimensional histogram implementation |
| `src/outputs/spherical_slice.cpp` | Fixed-radius binary angular-slice interpolation and writing |
