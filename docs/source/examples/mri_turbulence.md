# Example: MRI In A Shearing Box

The public tree contains three shearing-box magnetorotational-instability
(MRI) inputs using two generator-selection paths.

## Public Inputs

| Input deck | Dimension | Problem-generator path |
| --- | --- | --- |
| `inputs/shearing_box/mri3d_unstratified.athinput` | 3D | Built-in `pgen_name = mri3d` |
| `inputs/shearing_box/mri3d_stratified.athinput` | 3D | Built-in `pgen_name = mri3d`, with stratified setup controls in the deck |
| `inputs/shearing_box/mri2d.athinput` | 2D | Custom source `src/pgen/mri2d.cpp`; configure with `-DPROBLEM=mri2d` |

All three use `<shearing_box>` with `qshear` and `omega0`, `<mhd>`, and
x1 `shear_periodic` boundaries.

## Built-In 3D Validation Run

The unstratified 3D input can run with the default executable:

```bash
cmake -S . -B build
cmake --build build
./build/src/athena -i inputs/shearing_box/mri3d_unstratified.athinput \
  -d run-mri3d time/nlim=1
```

In a verified one-cycle run this deck wrote:

```text
run-mri3d/HGB.user.hst
run-mri3d/vtk/HGB.mhd_w_bcc.00000.vtk
run-mri3d/vtk/HGB.mhd_w_bcc.00001.vtk
```

The history stream uses `user_hist_only = true`, and
`src/pgen/tests/mri3d.cpp` enrolls the MRI history function.

## Custom 2D Status

The 2D input has no `pgen_name`; use its custom generator:

```bash
cmake -S . -B build-mri2d -DPROBLEM=mri2d
cmake --build build-mri2d
```

The shipped `inputs/shearing_box/mri2d.athinput` is **not currently a runnable
validation example**. In an audit run with this custom build, execution stops
at startup with:

```text
Shearing box source terms not enabled for mri2d problem
```

The input does define `<shearing_box>` and the MHD constructor creates its
shearing-box objects from that block. However, `src/pgen/mri2d.cpp` currently
checks `pmhd->psrc`, which is the generic `<mhd_srcterms>` object, rather than
the shearing-box object used by the built-in 3D generator. Until that
implementation/input mismatch is resolved and validated, use the verified
3D deck above for a public MRI run.

## Interpreting The Setup

The decks provide parameters such as `beta`, perturbation amplitude `amp`,
mode count `nwx`, and `ifield`; inspect the selected input and its generator
source together before changing them. The history data are intended for
monitoring MRI diagnostics, while VTK field dumps are used to examine the
spatial magnetic and fluid structure.

## See Also

- [Shearing Box](../modules/shearing_box.md)
- [MHD](../modules/mhd.md)
- [Problem Generators](../modules/pgen.md)
