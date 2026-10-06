# CGL Landau-fluid boundary contract

LF integration temporarily stores magnetic moment `mu = p_perp / |B|` in
`IAN`; ordinary active-CGL integration stores conservative anisotropy `A`.
This applies to explicit LF half-sweeps and RKL2 sweeps, on uniform and refined
meshes. `MHD::cgl_slot_representation` identifies the current representation.

## Fixed inflow

Populate `pbval_u->u_in` with the ordinary conserved CGL state, including total
energy and A. Populate `pbval_b->b_in` with the prescribed magnetic field.
The boundary implementation first fills face fields during an LF sweep, then
converts each inflow face's stored A to mu using the final cell-centered field.
Every face owns its conversion; overlapping corner slabs are not converted a
second time. Outflow, reflecting, and diode copies already carry the current
representation. Periodic communication also preserves that representation.

## User callbacks

A CGL user boundary must fill its face-field ghosts before converting primitive
pressures, and must supply cell-centered fields computed from those updated
faces to `CGLMHD::PrimToCons`. That method writes the representation currently
selected by the MHD module, including when its destination is a temporary
array. A callback may then copy the converted values into its owned ghost
cells. Calls to the scalar helper `SingleP2C_CGLMHD` still produce ordinary A;
code using that helper or directly writing conserved slots must explicitly
respect `cgl_slot_representation`.

Callbacks must be idempotent ghost fills. They must not advance physical state,
update active cells, or perform one-time side effects. On refined CGL-LF meshes,
physical boundaries are filled both before prolongation (for its stencil) and
after prolongation (to refresh corner donors before the next LF stencil).

`src/pgen/tests/cgl_lf_boundary.cpp` provides a complete `PrimToCons` callback
and an independently encoded fixed-inflow reference. The tests exercise explicit
LF and STS, fixed and user boundaries, 1D, mixed physical corners, primitive
SMR, and 3D, including CPU/GPU and one/four ranks. Uniform-state checks cover
active cells and the full one-cell LF halo. Nonuniform 1D checks compare
isotropic stored inflow with the prescribed-primitive callback.
