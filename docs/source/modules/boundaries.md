# Module: Boundary Values

Boundary-value code fills physical ghost zones and communicates data between
MeshBlocks. `Mesh` parses boundary names from `<mesh>`, while each active
physics module supplies its applicable physical-boundary implementation.

## Source Map

| Source | Responsibility |
| --- | --- |
| `src/bvals/bvals.*` | Boundary buffer classes and shared setup |
| `src/bvals/buffs_cc.cpp`, `buffs_fc.cpp` | Cell- and face-centered communication packing |
| `src/bvals/bvals_cc.cpp`, `bvals_fc.cpp` | Exchanges and unpacking |
| `src/bvals/flux_correct_cc.cpp`, `flux_correct_fc.cpp` | Multilevel flux correction |
| `src/bvals/prolongation.cpp`, `prolong_prims.cpp` | Multilevel prolongation |
| `src/bvals/physics/hydro_bcs.cpp` | Hydro physical boundaries |
| `src/bvals/physics/bfield_bcs.cpp` | Magnetic field physical boundaries |
| `src/bvals/physics/radiation_bcs.cpp` | Radiation physical boundaries |
| `src/bvals/physics/z4c_bcs.cpp`, `src/z4c/z4c_Sbc.cpp` | Z4c boundary handling |

## Parsed Boundary Names

`Mesh::GetBoundaryFlag()` recognizes:

| Input value | Role |
| --- | --- |
| `periodic` | Periodic pairing |
| `reflect` | Reflecting boundary |
| `inflow` | Inflow initialized through the selected problem setup |
| `outflow` | Outflow copy behavior |
| `diode` | Prevents inflow in fluid boundary implementations |
| `vacuum` | Vacuum behavior implemented for applicable state |
| `user` | User boundary hook |
| `shear_periodic` | Shearing-box x1 boundary |
| `undef` | Internal sentinel for inactive directions; not an active-face physical selection |

Not every physics module implements every parsed boundary type. In particular,
the current radiation physical-boundary source handles `outflow` and `inflow`
paths, while hydro and magnetic-field boundary sources include the wider fluid
set. Verify the active module's `src/bvals/physics/*_bcs.cpp` implementation
before selecting specialized boundaries.

```ini
<mesh>
ix1_bc = outflow
ox1_bc = outflow
ix2_bc = periodic
ox2_bc = periodic
```

## Constraints

- A periodic boundary on one face requires the opposite face in the same
  direction also to be periodic.
- `shear_periodic` is permitted only on both x1 faces, requires a
  `<shearing_box>` block, requires a multidimensional mesh, and cannot be
  combined with SMR/AMR in the public mesh constructor.
- `inflow` is not fully configured by a boundary string alone; the selected
  problem generator must initialize the inflow state.
- `user` requires a selected problem generator that enrolls the needed
  boundary callback.

## Multilevel And MHD Behavior

Cell-centered and face-centered communication use separate buffer paths. MHD
therefore communicates face-centered magnetic fields in addition to
cell-centered state and uses the face-centered flux-correction path on
multilevel meshes. The concrete prolongation and correction kernels are the
source of truth; no generic interpolation formula is promised by this page.

## See Also

- [Mesh](mesh.md)
- [Problem Generators](pgen.md)
- [Shearing Box](shearing_box.md)
