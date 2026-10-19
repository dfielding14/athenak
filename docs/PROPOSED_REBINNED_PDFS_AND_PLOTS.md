# Proposed rebinned PDFs and plot suite

Status: **implemented and release-validated**. The radius revision below passed
the 2026-06-09 golden/canary sequence, and the 69-product science campaign
completed all 161 inventory snapshots. Review changes made afterward require a
new release identity before any future production run.

## The short version

The current bins are not optimal, but the main problem is not simply that they
are too coarse. Several important axes have the wrong range:

- Spherical radius starts at 0.1 kpc, so it cannot show the resolved central
  region of a 4 pc or 8 pc simulation.
- Spherical radius stops at 190 kpc, so volume-weighted products put exactly
  55.11% of the cubic domain outside the plotted radial range.
- Small- and large-dust axes start at `1e-5`; their underflow contains a median
  56% of the mass.
- Vertical plots stop at 50 kpc in height and 30 kpc in cylindrical radius.
  Some vertical-inflow products exclude median fractions of 58% and 70% on
  those axes.
- The old 16-bin `|cos(theta)|` axis gives only about 2.1 bins within 30
  degrees of the z axis.
- Several legacy positive-velocity and north-hemisphere-only axes discard the
  opposite sign by construction.

The replacement therefore uses wider ranges, finer geometry bins, and
zero-aware symlog bins where zero is physically meaningful. It does **not**
blindly double every axis. That would waste memory on already adequate or
intrinsically narrow distributions.

Physical phase diagrams are never integrated over the whole domain in the
replacement. They retain coarse radius and `|cos(theta)|` axes, so every phase
diagram can be inspected by radial shell and angular sector. The replacement
also includes the requested four-dimensional temperature--radial-velocity
PDFs conditioned on radius and angle for the complete nonredundant family of
mass- and energy-flux weights.

The complete audit behind these statements is here:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/binning_audit/binning_audit.md
```

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/binning_audit/binning_audit.json
```

## Why the useless radius-costheta movie existed

`radius_costheta_small_dust_metallicity_mass_weighted.res_4pc_highmdot.mp4`
came from a generic pairwise marginal of the old four-dimensional `output35`
PDF. The plotter summed over small dust and metallicity and left radius versus
`|cos(theta)|`, weighted by total mass.

That is not a useful view of either dust or metallicity. It is a
mass-weighted geometry map. It will not exist in the replacement suite.

The replacement rule is:

> Never make a radius-versus-angle heatmap from a mass- or volume-weighted
> state PDF. Use that PDF only to condition a physical quantity on radius and
> angle. Radius-versus-angle heatmaps are allowed only for transport weights,
> where they show opening angle and flow geometry.

## Current and archived products

The earlier 76-product production trees were moved by same-filesystem directory
rename before this document was written. Their rebuilt PDFs, 24,633 files and
60 GiB, are retained at:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/delete_after_rebin_review_20260606/production
```

Their plots, 26,565 PNGs, 1,155 MP4s, and metadata are retained at:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild_plots/delete_after_rebin_review_20260606/production
```

The radius-revised 69-product science campaign subsequently populated:

```text
/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/production
```

A read-only 2026-07-10 audit found all 161 expected manifests, the exact same
69 product IDs in every snapshot, matching shard/output metadata, passing
geometry and volume guards, finite product sums, and nonempty headers and
payloads. The corresponding live production plot tree is absent; the archived
plot suite above still provides the existing visual baseline. Delete the
archived trees only after the new plots receive visual review.

## Global plot rules

1. No filename or visible label will use `outputXX`.
2. No mass- or volume-weighted radius-angle heatmap will be made.
3. Positive heatmaps will use white-low CMasher sequential maps such as
   `cmr.rainforest_r`, `cmr.chroma_r`, `cmr.torch_r`, and
   `cmr.voltage_r`.
4. Signed quantities will continue to use diverging maps.
5. Empty support remains light gray; small positive values fade toward white.
6. `R_cyl` is always the x axis and `|z|` is always the y axis.
7. Figure titles show fixed-width time relative to the first SN.
8. The six fixed radial shells remain 2, 4, 8, 16, 32, and 64 kpc.
9. The four angular selections remain all angles, polar 30 degrees,
   midplane 30 degrees, and the intermediate 30 degrees.
10. Movies are made only after the complete approved PNG suite passes audit.

## Proposed axis catalog

The aliases below define the exact bins used in every proposed PDF. Bin counts
exclude the underflow and overflow bins that the reducer adds automatically.

| Alias | Quantity | Scale | Interior bins | Range | Resolution or linthresh |
|---|---|---|---:|---|---|
| `R_FULL` | spherical radius | log | 361 | 0.064 to 262.144 kpc | 0.01001 dex/bin, factor 1.023 |
| `ABS_MU` | `|cos(theta)|` | linear | 64 | 0 to 1 | 0.015625/bin |
| `SIGNED_MU` | `cos(theta)` | linear | 128 | -1 to 1 | 0.015625/bin |
| `R_PHASE` | coarse spherical radius for 4D PDFs | log | 64 | 0.064 to 128 kpc | 0.05158 dex/bin, factor 1.126 |
| `ABS_MU_PHASE` | coarse `|cos(theta)|` for 4D PDFs | linear | 16 | 0 to 1 | 0.0625/bin |
| `Z_FULL` | `|z|` | log | 192 | 0.002 to 210 kpc | 0.02615 dex/bin, factor 1.062 |
| `RCYL_FULL` | cylindrical radius | log | 192 | 0.002 to 290 kpc | 0.02688 dex/bin, factor 1.064 |
| `P_KB` | pressure divided by `k_B` | log | 256 | `1e-6` to `1e10` K cm^-3 | 0.0625 dex/bin |
| `N_H` | hydrogen number density | log | 192 | `1e-10` to `1e6` cm^-3 | 0.08333 dex/bin |
| `T_K` | temperature | log | 192 | 1 to `1e9` K | 0.04688 dex/bin |
| `COOLING_RATE` | net cooling rate | symlog | 192 | `-1e-14` to `1e-14` erg s^-1 cm^3 | linthresh `1e-30` |
| `T_COOL` | net cooling time | symlog | 192 | `-1e8` to `1e8` Myr | linthresh 0.1 Myr |
| `V_R` | signed radial velocity | symlog | 192 | -10,000 to 10,000 km/s | linthresh 5 km/s |
| `V_TANGENTIAL` | signed polar or azimuthal velocity | symlog | 192 | -2,500 to 2,500 km/s | linthresh 5 km/s |
| `N_H_PHASE` | coarse phase-space density | log | 96 | `1e-10` to `1e6` cm^-3 | 0.16667 dex/bin |
| `T_PHASE` | coarse phase-space temperature | log | 96 | 1 to `1e9` K | 0.09375 dex/bin |
| `Z_GAS_PHASE` | coarse gas metallicity | positive symlog | 96 | 0 to 15 solar | linthresh `5e-7` solar |
| `DUST_PHASE` | coarse small- or large-grain dust ratio | positive symlog | 96 | 0 to 1 | linthresh `1e-8` |
| `T_TRANSPORT` | 4D transport temperature | log | 80 | 1 to `1e9` K | 0.1125 dex/bin |
| `V_R_TRANSPORT` | 4D signed radial velocity | symlog | 80 | -10,000 to 10,000 km/s | linthresh 5 km/s |
| `MACH` | total Mach number | positive symlog | 128 | 0 to 1,000 | linthresh `1e-3` |
| `RADIAL_MACH` | signed radial Mach number | symlog | 192 | -1,000 to 1,000 | linthresh `1e-2` |
| `ABS_RADIAL_MACH` | absolute radial Mach number | positive symlog | 128 | 0 to 1,000 | linthresh `1e-2` |
| `ABS_VR` | absolute radial velocity | positive symlog | 128 | 0 to 10,000 km/s | linthresh 1 km/s |
| `C_S` | sound speed | log | 128 | 0.1 to 10,000 km/s | 0.03906 dex/bin |
| `ENTROPY` | entropy proxy | log | 160 | `1e-8` to `1e9` code units | 0.10625 dex/bin |
| `J_Z` | specific angular momentum | symlog | 192 | `-1e6` to `1e6` kpc km/s | linthresh 10 kpc km/s |
| `Z_GAS` | gas metallicity | positive symlog | 96 | 0 to 15 solar | linthresh `5e-7` solar |
| `DUST` | small- or large-grain dust ratio | positive symlog | 128 | 0 to 1 | linthresh `1e-8` |
| `D_OVER_Z` | dust / total metal | linear | 96 | 0 to 1 | 0.01042/bin |
| `F_SMALL` | small-grain / total dust | linear | 96 | 0 to 1 | 0.01042/bin |

The positive-symlog axes are deliberate. They put exact zeros in an interior
bin and retain logarithmic resolution away from zero. This fixes the current
dust underflow without pretending that zero can be represented on a log axis.

## Plot recipes

The recipe codes keep the product catalog readable while defining exactly what
will be plotted from each PDF.

| Recipe | Plots made from the PDF |
|---|---|
| `S` | 2x2 quantity versus radius for the four angle selections; 3x2 quantity versus `|cos(theta)|` at the six fixed radial shells; 2x1 conditional mean/median profiles. **No radius-angle heatmap.** |
| `CP` | Four 3x2 physical-phase figures at the six fixed radial shells: one figure each for all angles, polar 30 degrees, midplane 30 degrees, and the intermediate band. No whole-domain phase diagram is made. |
| `O` | Radius-signed-angle transport heatmap; radius-temperature heatmap; 3x2 temperature-angle views at fixed radii; radial profiles by angle and thermal phase. |
| `G` | Radius-signed-angle transport heatmap; radial transport profiles for the four angle selections; angle profiles at the six fixed radial shells. |
| `TP` | Radius-temperature transport heatmap; radial profiles split by thermal phase; temperature distributions at the six fixed radial shells. |
| `TV` | Four 3x2 temperature--radial-velocity figures at the six fixed radial shells, one per angular selection; 2x2 temperature--radius and radial-velocity--radius figures for the four angular selections; folded radius--`|cos(theta)|` transport map. |
| `RM` | Radius-absolute-radial-Mach heatmap; sub/supersonic radial profiles; Mach distributions at fixed radial shells. |
| `VC` | Absolute-radial-velocity versus sound-speed heatmap with the Mach-one line; conditional median curves. |
| `VG` | `R_cyl` versus `|z|` heatmap; integrated vertical and cylindrical profiles. |
| `VP` | `|z|` versus temperature heatmap; height profiles split by thermal phase. |
| `J` | Radius-angular-momentum transport heatmap; radial angular-momentum transport profile. |

## Proposed PDF catalog

### Radial state

These PDFs answer how a physical state variable changes with radius and angle.
Their direct mass/volume radius-angle marginals are explicitly suppressed.

| PDF id | Weight | Axes | Dense bins | Plot recipe |
|---|---|---|---:|---|
| `science_r_theta_pressure_volume` | volume | `R_FULL x ABS_MU x P_KB` | 6,181,164 | `S` |
| `science_r_theta_density_volume` | volume | `R_FULL x ABS_MU x N_H` | 4,647,852 | `S` |
| `science_r_theta_temperature_volume` | volume | `R_FULL x ABS_MU x T_K` | 4,647,852 | `S` |
| `science_r_theta_temperature_edot_cool` | net cooling luminosity | `R_FULL x ABS_MU x T_K` | 4,647,852 | `S` |
| `science_r_theta_vr_volume` | volume | `R_FULL x ABS_MU x V_R` | 4,647,852 | `S` |
| `science_r_theta_vtheta_mass` | mass | `R_FULL x ABS_MU x V_TANGENTIAL` | 4,647,852 | `S` |
| `science_r_theta_vphi_mass` | mass | `R_FULL x ABS_MU x V_TANGENTIAL` | 4,647,852 | `S` |
| `science_r_theta_mach_mass` | mass | `R_FULL x ABS_MU x MACH` | 3,114,540 | `S` |
| `science_r_theta_radial_mach_mass` | mass | `R_FULL x ABS_MU x RADIAL_MACH` | 4,647,852 | `S` |
| `science_r_theta_entropy_mass` | mass | `R_FULL x ABS_MU x ENTROPY` | 3,881,196 | `S` |
| `science_r_theta_jz_mass` | mass | `R_FULL x ABS_MU x J_Z` | 4,647,852 | `S` |
| `science_r_theta_metallicity_mass` | mass | `R_FULL x ABS_MU x Z_GAS` | 2,347,884 | `S` |
| `science_r_theta_small_dust_mass` | mass | `R_FULL x ABS_MU x DUST` | 3,114,540 | `S` |
| `science_r_theta_large_dust_mass` | mass | `R_FULL x ABS_MU x DUST` | 3,114,540 | `S` |
| `science_r_theta_edot_cool_volume` | volume | `R_FULL x ABS_MU x COOLING_RATE` | 4,647,852 | `S` |
| `science_r_theta_tcool_volume` | volume | `R_FULL x ABS_MU x T_COOL` | 4,647,852 | `S` |

### Physical phase structure

The previous draft was wrong here: its phase PDFs had no radial or angular
conditioning. The replacement uses four-dimensional PDFs of the form
`R_PHASE x ABS_MU_PHASE x Q1 x Q2`. The reducer supports at most four axes, so
retaining radius and angle means each PDF must contain one physical pair
rather than three physical variables.

Every product below has 64 radial bins and 16 angular bins. The polar
30-degree cap spans about 2.14 coarse angular bins. That is intentionally
coarser than the radial-state products, but it preserves the dependence that a
whole-domain phase diagram would erase.

| PDF id | Weight | Axes | Dense bins | Plot recipe |
|---|---|---|---:|---|
| `science_phase_density_temperature_volume` | volume | `R_PHASE x ABS_MU_PHASE x N_H_PHASE x T_PHASE` | 11,409,552 | `CP` |
| `science_phase_density_temperature_mass` | mass | `R_PHASE x ABS_MU_PHASE x N_H_PHASE x T_PHASE` | 11,409,552 | `CP` |
| `science_phase_temperature_metallicity_mass` | mass | `R_PHASE x ABS_MU_PHASE x T_PHASE x Z_GAS_PHASE` | 11,409,552 | `CP` |
| `science_phase_temperature_small_dust_mass` | mass | `R_PHASE x ABS_MU_PHASE x T_PHASE x DUST_PHASE` | 11,409,552 | `CP` |
| `science_phase_temperature_large_dust_mass` | mass | `R_PHASE x ABS_MU_PHASE x T_PHASE x DUST_PHASE` | 11,409,552 | `CP` |
| `science_phase_metallicity_small_dust_mass` | mass | `R_PHASE x ABS_MU_PHASE x Z_GAS_PHASE x DUST_PHASE` | 11,409,552 | `CP` |
| `science_phase_metallicity_large_dust_mass` | mass | `R_PHASE x ABS_MU_PHASE x Z_GAS_PHASE x DUST_PHASE` | 11,409,552 | `CP` |
| `science_phase_dust_survival` | total metal mass | `R_PHASE x ABS_MU_PHASE x T_PHASE x D_OVER_Z` | 11,409,552 | `CP` |
| `science_phase_grain_processing` | dust mass | `R_PHASE x ABS_MU_PHASE x T_PHASE x F_SMALL` | 11,409,552 | `CP` |
| `science_phase_density_temperature_edot_cool` | net cooling luminosity | `R_PHASE x ABS_MU_PHASE x N_H_PHASE x T_PHASE` | 11,409,552 | `CP` |

### Opening angle and radial transport geometry

Radius-angle heatmaps are retained here because the weights are transport
rates, so the maps show genuine flow geometry.

The two 3D outflow PDFs also supply the outflow radius-angle and
radius-temperature marginals. Separate duplicate 2D outflow PDFs are not
needed.

| PDF id | Weight | Axes | Dense bins | Plot recipe |
|---|---|---|---:|---|
| `science_phase_opening_angle_mdot_out` | outflow mass flux | `R_FULL x SIGNED_MU x T_K` | 9,154,860 | `O` |
| `science_phase_opening_angle_edot_out` | outflow energy flux | `R_FULL x SIGNED_MU x T_K` | 9,154,860 | `O` |
| `science_transport_geometry_mdot_in_abs` | inflow mass-flux magnitude | `R_FULL x SIGNED_MU` | 47,190 | `G` |
| `science_transport_geometry_edot_in_abs` | inflow energy-flux magnitude | `R_FULL x SIGNED_MU` | 47,190 | `G` |
| `science_transport_geometry_edot_kin_out` | outflow kinetic-energy flux | `R_FULL x SIGNED_MU` | 47,190 | `G` |
| `science_transport_geometry_edot_kin_in_abs` | inflow kinetic-energy-flux magnitude | `R_FULL x SIGNED_MU` | 47,190 | `G` |
| `science_transport_geometry_edot_th_out` | outflow enthalpy flux | `R_FULL x SIGNED_MU` | 47,190 | `G` |
| `science_transport_geometry_edot_th_in_abs` | inflow enthalpy-flux magnitude | `R_FULL x SIGNED_MU` | 47,190 | `G` |
| `science_transport_geometry_ram_out` | outflow ram flux | `R_FULL x SIGNED_MU` | 47,190 | `G` |

### Four-dimensional thermokinematic transport

These are the requested temperature--radial-velocity--radius--angle PDFs.
They use folded `|cos(theta)|`; the signed-angle products above remain
necessary for north--south asymmetry.

The inflow products store positive magnitudes. Separate negative-valued
`mdot_in` and `edot_in` products would contain the same information with the
opposite sign and are therefore not duplicated.

| PDF id | Weight | Axes | Dense bins | Plot recipe |
|---|---|---|---:|---|
| `science_transport_thermokinematic_mdot_net` | signed net mass flux | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_mdot_out` | outflow mass flux | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_mdot_in_abs` | inflow mass-flux magnitude | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_net` | signed net total-energy flux | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_out` | outflow total-energy flux | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_in_abs` | inflow total-energy-flux magnitude | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_kin_net` | signed net kinetic-energy flux | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_kin_out` | outflow kinetic-energy flux | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_kin_in_abs` | inflow kinetic-energy-flux magnitude | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_th_net` | signed net enthalpy flux | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_th_out` | outflow enthalpy flux | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_th_in_abs` | inflow enthalpy-flux magnitude | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |
| `science_transport_thermokinematic_edot_cool` | net cooling luminosity | `R_PHASE x ABS_MU_PHASE x T_TRANSPORT x V_R_TRANSPORT` | 7,988,112 | `TV` |

### Composition transport by thermal phase

The new 4D thermokinematic PDFs supply all gas mass- and energy-flux
temperature marginals. The products below retain the independent metal and
dust weights.

| PDF id | Weight | Axes | Dense bins | Plot recipe |
|---|---|---|---:|---|
| `science_gas_metal_mdot_out` | outflow gas-metal flux | `R_FULL x T_K` | 70,422 | `TP` |
| `science_gas_metal_mdot_in_abs` | inflow gas-metal-flux magnitude | `R_FULL x T_K` | 70,422 | `TP` |
| `science_total_metal_mdot_out` | outflow total-metal flux | `R_FULL x T_K` | 70,422 | `TP` |
| `science_total_metal_mdot_in_abs` | inflow total-metal-flux magnitude | `R_FULL x T_K` | 70,422 | `TP` |
| `science_dust_mdot_out` | outflow dust flux | `R_FULL x T_K` | 70,422 | `TP` |
| `science_dust_mdot_in_abs` | inflow dust-flux magnitude | `R_FULL x T_K` | 70,422 | `TP` |

### Radial dynamics

| PDF id | Weight | Axes | Dense bins | Plot recipe |
|---|---|---|---:|---|
| `science_radial_mach_mass` | mass | `R_FULL x ABS_RADIAL_MACH` | 47,190 | `RM` |
| `science_radial_mach_mdot_out` | outflow mass flux | `R_FULL x ABS_RADIAL_MACH` | 47,190 | `RM` |
| `science_radial_mach_mdot_in_abs` | inflow mass-flux magnitude | `R_FULL x ABS_RADIAL_MACH` | 47,190 | `RM` |
| `science_vr_sound_speed_mdot_out` | outflow mass flux | `ABS_VR x C_S` | 16,900 | `VC` |
| `science_vr_sound_speed_mdot_in_abs` | inflow mass-flux magnitude | `ABS_VR x C_S` | 16,900 | `VC` |
| `science_vr_sound_speed_edot_out` | outflow energy flux | `ABS_VR x C_S` | 16,900 | `VC` |
| `science_vr_sound_speed_ram_out` | outflow ram flux | `ABS_VR x C_S` | 16,900 | `VC` |

### Vertical transport

| PDF id | Weight | Axes | Dense bins | Plot recipe |
|---|---|---|---:|---|
| `science_vertical_geometry_mdot_out` | vertical outflow mass flux | `RCYL_FULL x Z_FULL` | 37,636 | `VG` |
| `science_vertical_geometry_mdot_in_abs` | vertical inflow mass-flux magnitude | `RCYL_FULL x Z_FULL` | 37,636 | `VG` |
| `science_vertical_geometry_edot_out` | vertical outflow energy flux | `RCYL_FULL x Z_FULL` | 37,636 | `VG` |
| `science_vertical_geometry_ram_out` | vertical outflow ram flux | `RCYL_FULL x Z_FULL` | 37,636 | `VG` |
| `science_vertical_phase_mdot_out` | vertical outflow mass flux | `Z_FULL x T_K` | 37,636 | `VP` |
| `science_vertical_phase_edot_out` | vertical outflow energy flux | `Z_FULL x T_K` | 37,636 | `VP` |

### Angular-momentum transport

The separate mass-weighted radius-angular-momentum PDF is removed because it is
already a marginal of `science_r_theta_jz_mass`.

| PDF id | Weight | Axes | Dense bins | Plot recipe |
|---|---|---|---:|---|
| `science_angular_momentum_mdot_out` | outflow mass flux | `R_FULL x J_Z` | 70,422 | `J` |
| `science_angular_momentum_mdot_in_abs` | inflow mass-flux magnitude | `R_FULL x J_Z` | 70,422 | `J` |

## Products explicitly removed

The replacement catalog contains 69 descriptive science PDFs. These current
products are intentionally not recalculated:

- All legacy `output23` through `output36` products.
- Redundant 2D physical marginals:
  `science_density_temperature_volume`,
  `science_temperature_metallicity_mass`,
  `science_temperature_small_dust_mass`, and
  `science_temperature_large_dust_mass`.
- Folded radius-angle transport products:
  `science_r_theta_mdot_out`, `science_r_theta_mdot_in`,
  `science_r_theta_edot_out`, `science_r_theta_edot_in`,
  `science_r_theta_edot_kin`, and `science_r_theta_edot_th`.
- Duplicate outflow marginals already supplied by the 3D opening-angle PDFs:
  `science_transport_geometry_mdot_out`,
  `science_transport_geometry_edot_out`,
  `science_transport_phase_mdot_out`, and
  `science_transport_phase_edot_out`.
- Low-dimensional gas transport-phase products now supplied by the 4D
  thermokinematic PDFs:
  `science_transport_phase_mdot_in_abs`,
  `science_transport_phase_edot_in_abs`,
  `science_transport_phase_edot_kin_out`,
  `science_transport_phase_edot_kin_in_abs`,
  `science_transport_phase_edot_th_out`, and
  `science_transport_phase_edot_th_in_abs`.
- Global physical-phase products and the old radius-only density--temperature
  products are replaced by the conditioned `R_PHASE x ABS_MU_PHASE` suite.
- Redundant marginals `science_entropy_mass` and
  `science_angular_momentum_mass`.

## What the audit supports

This is the part I trust:

- `R_FULL` fixes both the missing central region and the outer-domain loss.
- `ABS_MU=64` and `SIGNED_MU=128` give the same angular width and resolve the
  polar 30-degree region with about 8.6 bins instead of 2.1.
- `Z_FULL` and `RCYL_FULL` cover the full simulation domain instead of
  truncating vertical inflow and late-time outflow.
- Zero-aware dust and metallicity axes fix a real underflow problem.
- Physical temperature and density axes have better ranges than the old code
  variable axes.
- Every physical phase diagram now retains radius and angular sector.
- The 4D transport suite retains temperature, signed radial velocity, radius,
  and angle simultaneously for all nonredundant mass- and energy-flux
  variants.

The completed golden and production runs resolved the release-level questions,
but these interpretation checks remain for plot review:

- Quantify exact-zero versus merely tiny dust values from the new
  positive-symlog bins; the old underflow bin could not separate them.
- Use the rebuilt metadata to measure the remaining underflow and overflow in
  each scientific view; passing the reducer checks does not decide whether a
  presentation range is optimal.
- Every subsequent source or product-catalog change creates a new executable
  identity and must repeat the complete golden, validator, and canary sequence.

## Approved implementation decisions

1. The 69-PDF catalog is the current target.
2. `R_FULL=361` over 0.064 to 262.144 kpc and `R_PHASE=64` over 0.064 to 128
   kpc are the current target bins.
3. Only transport-weighted products may make radius-angle heatmaps.
4. Duplicate outflow 2D views are derived from the two 3D opening-angle PDFs
   instead of stored separately.
5. The radius-revised release gate is 6 GiB/rank, because the requested
   `R_FULL` grid pushes the planned peak buffer to 4.806 GiB/rank.

## Compute and storage cost

The pre-run estimate used the earlier completed 161-snapshot campaign as its
empirical baseline. That campaign processed 11.100 trillion cells, read 355.2
TB of source payload, and consumed 67.54 actual reducer node-hours. The
implemented 69-product campaign processed the same cell and byte totals and
used 125.55 reducer node-hours, within the predicted 100--160 node-hour range.

| Quantity | Historical 76-product campaign | Implemented rebinned campaign |
|---|---:|---:|
| PDFs per snapshot | 76 | 69 |
| Atomic histogram updates per cell | 76 | 69 |
| Dense bins per rank | 43,593,352 | 305,811,772 |
| Histogram memory per rank | 0.325 GiB | 2.278 GiB |
| Planned peak buffer per rank | 1.150 GiB with `CHUNK_BLOCKS=32` | 4.806 GiB with `CHUNK_BLOCKS=16` |
| Declared release memory budget | 1.5 GiB/rank | 6.0 GiB/rank |
| Dense PDF payload for all 161 snapshots | 52.3 GiB theoretical, 60 GiB tree | 366.8 GiB theoretical, about 410 GiB tree |
| Plot groups per snapshot | 165 | about 315 |
| Full PNG count | 26,565 | about 50,715 |
| Movie count | 1,155 | about 2,205 |

The new run performs 9% fewer per-cell product updates, and it still reads
the 355.2 TB source only once. Its dense MPI reductions and writes are about
4.66 times larger, however. The old 1.5 GiB/rank release gate cannot be reused:
the new catalog needs a predeclared 6.0 GiB/rank gate and must prove the actual
peak and runtime on the new golden before production submission.

The previous golden timings show that input scanning and per-cell updates
dominate, while the larger reductions and output add real but secondary cost.
The pre-run estimate was:

- Expected full production reducer cost: **100 to 160 Frontier node-hours**;
  actual: **125.55 Frontier node-hours**.
- Conservative reservation for production: **200 Frontier node-hours**.
- Repeated goldens and canaries: **about 2 to 5 Frontier node-hours**.
- Plot and movie campaign: **about 1 to 2 service-node-hours**.
- New durable storage before deleting the archive: **about 440 GiB** for
  rebinned PDFs, plots, movies, and metadata.

No dollar estimate is given because Frontier allocation accounting is in
node-hours rather than a project-specific public dollar rate.
