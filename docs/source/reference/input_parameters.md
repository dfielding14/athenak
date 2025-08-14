# AthenaK Input Parameters Reference
Complete list of all input parameters by block, extracted from source code.

## Input Block: `<coord>`
**Used by**: coordinates.cpp, units.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `a` | Real | required | coordinates.cpp:L47 |
| `dexcise` | Real | required | coordinates.cpp:L57 |
| `excise` | bool | true | coordinates.cpp:L48 |
| `excise_lapse` | Real | 0.25 | coordinates.cpp:L72 |
| `excision_scheme` | string | fixed | coordinates.cpp:L67 |
| `flux_excise_r` | Real | 1.0 | coordinates.cpp:L61 |
| `general_rel` | bool | false | units.cpp:L24 |
| `minkowski` | bool | false | coordinates.cpp:L45 |
| `pexcise` | Real | required | coordinates.cpp:L58 |
| `special_rel` | bool | false | coordinates.cpp:L32 |

## Input Block: `<eos>`
**Used by**: dynamo.cpp, sfb_turb.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `gamma` | Real | 5.0/3.0 | dynamo.cpp:L55 |
| `iso_sound_speed` | Real | 1.0 | sfb_turb.cpp:L50 |

## Input Block: `<hydro>`
**Used by**: hydro.cpp, ideal_grhyd.cpp, isothermal_hyd.cpp, rt.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `const_accel_val` | Real | required | rt.cpp:L85 |
| `eos` | string | required | hydro.cpp:L42 |
| `fofc` | bool | false | hydro.cpp:L134 |
| `gamma` | Real | required | ideal_grhyd.cpp:L28 |
| `gamma_max` | Real | (FLT_MAX | ideal_grhyd.cpp:L32 |
| `iso_sound_speed` | Real | required | isothermal_hyd.cpp:L21 |
| `nscalars` | int | 0 | hydro.cpp:L72 |
| `reconstruct` | string | plm | hydro.cpp:L137 |
| `rsolver` | string | required | hydro.cpp:L183 |

## Input Block: `<ion-neutral>`
**Used by**: ion-neutral.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `drag_coeff` | Real | required | ion-neutral.cpp:L28 |
| `ionization_coeff` | Real | 0.0 | ion-neutral.cpp:L29 |
| `recombination_coeff` | Real | 0.0 | ion-neutral.cpp:L30 |

## Input Block: `<job>`
**Used by**: outputs.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `basename` | string | required | outputs.cpp:L93 |

## Input Block: `<mesh>`
**Used by**: mesh.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `ix1_bc` | string | required | mesh.cpp:L80 |
| `ix2_bc` | string | required | mesh.cpp:L119 |
| `ix3_bc` | string | required | mesh.cpp:L146 |
| `nghost` | Real | 2 | mesh.cpp:L63 |
| `nx1` | int | required | mesh.cpp:L64 |
| `nx2` | int | required | mesh.cpp:L65 |
| `nx3` | int | required | mesh.cpp:L66 |
| `ox1_bc` | string | required | mesh.cpp:L81 |
| `ox2_bc` | string | required | mesh.cpp:L120 |
| `ox3_bc` | string | required | mesh.cpp:L147 |
| `x1max` | Real | required | mesh.cpp:L57 |
| `x1min` | Real | required | mesh.cpp:L56 |
| `x2max` | Real | required | mesh.cpp:L59 |
| `x2min` | Real | required | mesh.cpp:L58 |
| `x3max` | Real | required | mesh.cpp:L61 |
| `x3min` | Real | required | mesh.cpp:L60 |

## Input Block: `<mesh_refinement>`
**Used by**: build_tree.cpp, mesh.cpp, mesh_refinement.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `ddens_max` | Real | required | mesh_refinement.cpp:L69 |
| `dens_max` | Real | required | mesh_refinement.cpp:L65 |
| `dpres_max` | Real | required | mesh_refinement.cpp:L73 |
| `dvel_max` | Real | required | mesh_refinement.cpp:L77 |
| `max_nmb_per_rank` | Real | required | build_tree.cpp:L479 |
| `ncycle_check` | Real | 1 | mesh_refinement.cpp:L57 |
| `num_levels` | int | 1 | build_tree.cpp:L383 |
| `prolong_primitives` | bool | required | mesh_refinement.cpp:L61 |
| `refinement` | string | required | mesh.cpp:L175 |
| `refinement_interval` | Real | 5 | mesh_refinement.cpp:L58 |

## Input Block: `<meshblock>`
**Used by**: mesh.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `nx1` | int | mesh_indcs.nx1 | mesh.cpp:L245 |
| `nx2` | int | mesh_indcs.nx2 | mesh.cpp:L247 |
| `nx3` | int | mesh_indcs.nx3 | mesh.cpp:L252 |

## Input Block: `<mhd>`
**Used by**: dyn_grmhd.cpp, dyngr_tov.cpp, ideal_srmhd.cpp, isothermal_mhd.cpp, mhd.cpp, resistivity.cpp, rt.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `const_accel_val` | Real | required | rt.cpp:L91 |
| `dmp_M` | Real | 1.2 | dyn_grmhd.cpp:L122 |
| `dyn_eos` | string | required | dyn_grmhd.cpp:L57 |
| `dyn_error` | string | required | dyn_grmhd.cpp:L58 |
| `dyn_scratch` | int | 0 | dyn_grmhd.cpp:L120 |
| `enforce_maximum` | bool | true | dyn_grmhd.cpp:L121 |
| `eos` | string | required | mhd.cpp:L60 |
| `fixed` | bool | false | dyn_grmhd.cpp:L124 |
| `fofc` | bool | false | mhd.cpp:L176 |
| `fofc_method` | string | llf | dyn_grmhd.cpp:L109 |
| `gamma` | Real | required | dyngr_tov.cpp:L114 |
| `gamma_max` | Real | (FLT_MAX | ideal_srmhd.cpp:L26 |
| `iso_sound_speed` | Real | required | isothermal_mhd.cpp:L21 |
| `nscalars` | int | 0 | mhd.cpp:L92 |
| `ohmic_resistivity` | Real | required | resistivity.cpp:L26 |
| `reconstruct` | string | plm | mhd.cpp:L179 |
| `rsolver` | string | required | dyn_grmhd.cpp:L98 |

## Input Block: `<particles>`
**Used by**: outputs.cpp, part_static_turb.cpp, particles.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `assign_tag` | string | index_order | particles.cpp:L146 |
| `cfl_part` | Real | 0.05 | part_static_turb.cpp:L53 |
| `interpolation` | string | tsc | particles.cpp:L83 |
| `mass_log_spacing` | Real | 1.0 | part_static_turb.cpp:L57 |
| `min_mass` | Real | 1.0 | part_static_turb.cpp:L56 |
| `n_outputs` | int | required | outputs.cpp:L66 |
| `nspecies` | int | 1 | particles.cpp:L37 |
| `particle_type` | string | required | particles.cpp:L69 |
| `ppc` | Real | 1.0 | particles.cpp:L36 |
| `pusher` | string | required | particles.cpp:L82 |

## Input Block: `<potential>`
**Used by**: cgm_cooling_flow.cpp, cgm_cooling_flow_magnetized.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `mass_gal` | Real | required | cgm_cooling_flow_magnetized.cpp:L102 |
| `r_200` | Real | required | cgm_cooling_flow_magnetized.cpp:L104 |
| `r_scale` | Real | required | cgm_cooling_flow_magnetized.cpp:L100 |
| `rho_mean` | Real | required | cgm_cooling_flow_magnetized.cpp:L105 |
| `rho_scale` | Real | required | cgm_cooling_flow_magnetized.cpp:L101 |
| `scale_gal` | Real | required | cgm_cooling_flow_magnetized.cpp:L103 |
| `z_gal` | Real | required | cgm_cooling_flow.cpp:L91 |

## Input Block: `<problem>`
**Used by**: blast.cpp, cgm_cooling_flow_amr.cpp, cgm_cooling_flow_magnetized.cpp, cshock.cpp, current_sheet.cpp, dynamo.cpp, dyngr_tov.cpp, field_loop.cpp, gr_torus.cpp, hydrostatic_1d.cpp, hydrostatic_3d.cpp, kh.cpp, lorene_bns.cpp, mri2d.cpp, part_static_turb.cpp, pgen.cpp, rad_beam.cpp, rad_diffusion.cpp, rad_relax.cpp, rad_snake.cpp, rst_prtcl.cpp, rt.cpp, rt_mhd_amr.cpp, sfb_turb.cpp, sgrid_bns.cpp, shock_cloud.cpp, shwave.cpp, slotted_cyl.cpp, turb.cpp, turb_mhd_amr_wave.cpp, z4c_one_puncture.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `180rotation` | int | 0 | sgrid_bns.cpp:L163 |
| `A_snake` | Real | 0.1 | blast.cpp:L161 |
| `B0z` | Real | 1.0 | part_static_turb.cpp:L54 |
| `Interpolate_make_finer_grid2` | bool | false | sgrid_bns.cpp:L493 |
| `Interpolate_max_xyz_diff` | Real | 0.0 | sgrid_bns.cpp:L495 |
| `Interpolate_verbose` | bool | false | sgrid_bns.cpp:L491 |
| `Mach` | Real | required | shock_cloud.cpp:L38 |
| `Omega` | Real | required | sgrid_bns.cpp:L155 |
| `a0` | Real | 1.0 | current_sheet.cpp:L42 |
| `a_char` | Real | 0.01 | kh.cpp:L37 |
| `a_warp` | Real | 6.0 | blast.cpp:L160 |
| `amp` | Real | 0.01 | rt_mhd_amr.cpp:L44 |
| `b0` | Real | 0.1 | rt_mhd_amr.cpp:L46 |
| `b_amb` | Real | 0.1 | blast.cpp:L140 |
| `b_norm` | Real | 0.0 | dyngr_tov.cpp:L506 |
| `beta` | Real | 1.0 | sfb_turb.cpp:L51 |
| `bg` | Real | 0. | current_sheet.cpp:L43 |
| `by0` | Real | required | cshock.cpp:L80 |
| `center_x1` | Real | 0.50 | slotted_cyl.cpp:L76 |
| `center_x2` | Real | 0.75 | slotted_cyl.cpp:L77 |
| `chakrabarti_torus` | bool | false | gr_torus.cpp:L236 |
| `coordinates` | string | cartesian | blast.cpp:L141 |
| `d0` | Real | required | shwave.cpp:L55 |
| `d_i` | Real | 1.0 | turb.cpp:L117 |
| `d_n` | Real | 1.0 | turb.cpp:L118 |
| `datadir` | string | required | sgrid_bns.cpp:L624 |
| `di0` | Real | required | cshock.cpp:L74 |
| `di_amb` | Real | 1.0 | blast.cpp:L136 |
| `dir_1` | Real | required | rad_beam.cpp:L56 |
| `dir_2` | Real | required | rad_beam.cpp:L57 |
| `dir_3` | Real | required | rad_beam.cpp:L58 |
| `disk_profile_file` | string | required | cgm_cooling_flow_amr.cpp:L117 |
| `dn0` | Real | required | cshock.cpp:L75 |
| `dn_amb` | Real | 1.0 | blast.cpp:L133 |
| `dr` | Real | required | dyngr_tov.cpp:L695 |
| `drat` | Real | 1.0 | blast.cpp:L139 |
| `drho_rho0` | Real | 0.0 | kh.cpp:L43 |
| `ecc` | Real | required | sgrid_bns.cpp:L156 |
| `epsb` | Real | 0.0 | current_sheet.cpp:L45 |
| `epsv` | Real | 0.0 | current_sheet.cpp:L46 |
| `erad` | Real | required | rad_relax.cpp:L46 |
| `error_output` | int | required | shwave.cpp:L80 |
| `field_config` | int | 0 | turb_mhd_amr_wave.cpp:L63 |
| `fm_torus` | bool | false | gr_torus.cpp:L235 |
| `grav_acc` | Real | 1.0 | hydrostatic_1d.cpp:L40 |
| `ifield` | int | 1 | mri2d.cpp:L54 |
| `initial_data_file` | string | required | lorene_bns.cpp:L83 |
| `inner_radius` | Real | required | blast.cpp:L130 |
| `ipert` | int | 1 | rt_mhd_amr.cpp:L47 |
| `iprob` | int | required | rt.cpp:L68 |
| `isotropic` | bool | false | dyngr_tov.cpp:L708 |
| `k_snake` | Real | 2.0 | blast.cpp:L162 |
| `kappa` | Real | required | dyngr_tov.cpp:L113 |
| `keep_sgrid_output` | bool | false | sgrid_bns.cpp:L489 |
| `kval` | Real | 1.0 | current_sheet.cpp:L47 |
| `magindex` | Real | 2 | dyngr_tov.cpp:L508 |
| `mass_loss_radius` | Real | required | cgm_cooling_flow_magnetized.cpp:L108 |
| `mass_loss_rate` | Real | required | cgm_cooling_flow_magnetized.cpp:L107 |
| `minkowski` | bool | false | dyngr_tov.cpp:L356 |
| `n_param` | Real | 0.0 | gr_torus.cpp:L233 |
| `ng` | Real | 1.0 | current_sheet.cpp:L40 |
| `nhigh_ICs` | int | 4 | dynamo.cpp:L98 |
| `nlow_ICs` | int | 1 | dynamo.cpp:L97 |
| `npoints` | Real | required | dyngr_tov.cpp:L694 |
| `nu` | Real | required | rad_diffusion.cpp:L61 |
| `nwx` | int | required | shwave.cpp:L76 |
| `nwy` | int | required | shwave.cpp:L77 |
| `nwz` | int | required | shwave.cpp:L78 |
| `omega` | Real | 1.0 | slotted_cyl.cpp:L79 |
| `omega_x1` | Real | 0.50 | slotted_cyl.cpp:L80 |
| `omega_x2` | Real | 0.50 | slotted_cyl.cpp:L81 |
| `outdir` | string | SGRID | sgrid_bns.cpp:L496 |
| `outer_radius` | Real | required | blast.cpp:L129 |
| `p0` | Real | 1.0 | sfb_turb.cpp:L53 |
| `p_pert` | Real | 0.0 | dyngr_tov.cpp:L707 |
| `pcut` | Real | 1e-6 | dyngr_tov.cpp:L507 |
| `pert` | Real | 1.0e-4 | cshock.cpp:L81 |
| `pert_amp` | Real | 0.0 | gr_torus.cpp:L239 |
| `pgas0` | Real | 1.0 | turb_mhd_amr_wave.cpp:L55 |
| `pgas_min` | Real | required | gr_torus.cpp:L225 |
| `pgas_pow` | Real | required | gr_torus.cpp:L226 |
| `pgen_name` | string | none | pgen.cpp:L620 |
| `pi_amb` | Real | 1.0 | blast.cpp:L135 |
| `pn_amb` | Real | 1.0 | blast.cpp:L132 |
| `pos_1` | Real | required | rad_snake.cpp:L190 |
| `pos_2` | Real | required | rad_snake.cpp:L191 |
| `pos_3` | Real | required | rad_snake.cpp:L192 |
| `potential_beta_min` | Real | 100.0 | gr_torus.cpp:L513 |
| `potential_cutoff` | Real | 0.2 | gr_torus.cpp:L514 |
| `potential_falloff` | Real | 0.0 | gr_torus.cpp:L526 |
| `potential_r_pow` | Real | 0.0 | gr_torus.cpp:L527 |
| `potential_rho_pow` | Real | 1.0 | gr_torus.cpp:L528 |
| `prat` | Real | required | blast.cpp:L138 |
| `press` | Real | 1.0 | kh.cpp:L42 |
| `profile_file` | string | required | cgm_cooling_flow_magnetized.cpp:L114 |
| `prograde` | bool | true | gr_torus.cpp:L234 |
| `prtcl_res_file` | string | required | part_static_turb.cpp:L77 |
| `prtcl_rst_flag` | int | 0 | rst_prtcl.cpp:L42 |
| `punc_ADM_mass` | Real | 1. | z4c_one_puncture.cpp:L87 |
| `punc_center_x1` | Real | 0. | z4c_one_puncture.cpp:L88 |
| `punc_center_x2` | Real | 0. | z4c_one_puncture.cpp:L89 |
| `punc_center_x3` | Real | 0. | z4c_one_puncture.cpp:L90 |
| `r_circ` | Real | required | cgm_cooling_flow_magnetized.cpp:L110 |
| `r_core` | Real | 0.1 | sfb_turb.cpp:L56 |
| `r_edge` | Real | required | gr_torus.cpp:L231 |
| `r_outer` | Real | 1.0 | sfb_turb.cpp:L57 |
| `r_peak` | Real | required | gr_torus.cpp:L232 |
| `r_scale` | Real | 1.0 | hydrostatic_3d.cpp:L40 |
| `rad` | Real | 0.0 | field_loop.cpp:L57 |
| `radius` | Real | 0.15 | slotted_cyl.cpp:L75 |
| `rdot` | Real | required | sgrid_bns.cpp:L160 |
| `rho0` | Real | 1.0 | sfb_turb.cpp:L52 |
| `rho1` | Real | 1.0 | kh.cpp:L39 |
| `rho_max` | Real | required | gr_torus.cpp:L230 |
| `rho_min` | Real | required | gr_torus.cpp:L223 |
| `rho_outer` | Real | 0.1 | sfb_turb.cpp:L58 |
| `rho_pow` | Real | required | gr_torus.cpp:L224 |
| `rhoc` | Real | required | dyngr_tov.cpp:L693 |
| `s_height` | Real | 0.25 | slotted_cyl.cpp:L84 |
| `s_width` | Real | 0.05 | slotted_cyl.cpp:L83 |
| `sigma` | Real | required | kh.cpp:L35 |
| `smooth_interface` | bool | false | rt.cpp:L70 |
| `snake_kym` | Real | required | rad_snake.cpp:L58 |
| `snake_mag` | Real | required | rad_snake.cpp:L57 |
| `snake_tet` | bool | false | rad_snake.cpp:L59 |
| `spread` | Real | required | rad_snake.cpp:L194 |
| `table` | string | required | dyngr_tov.cpp:L188 |
| `tangled_ICs` | bool | false | dynamo.cpp:L56 |
| `temp` | Real | required | rad_relax.cpp:L47 |
| `tilt_angle` | Real | 0.0 | gr_torus.cpp:L227 |
| `use_pcut_rel` | bool | false | dyngr_tov.cpp:L512 |
| `user_hist` | bool | false | pgen.cpp:L145 |
| `user_srcs` | bool | false | pgen.cpp:L144 |
| `v1` | Real | required | rad_diffusion.cpp:L55 |
| `v_circ` | Real | required | cgm_cooling_flow_magnetized.cpp:L111 |
| `v_pert` | Real | 0.0 | dyngr_tov.cpp:L706 |
| `verbose` | bool | 0 | sgrid_bns.cpp:L611 |
| `vertical_field` | bool | false | gr_torus.cpp:L516 |
| `vix0` | Real | required | cshock.cpp:L76 |
| `viy0` | Real | required | cshock.cpp:L78 |
| `vnx0` | Real | required | cshock.cpp:L77 |
| `vny0` | Real | required | cshock.cpp:L79 |
| `vshear` | Real | required | kh.cpp:L36 |
| `vx` | Real | 1./1.2 | field_loop.cpp:L328 |
| `vx0` | Real | 0.0 | turb_mhd_amr_wave.cpp:L56 |
| `vy` | Real | 1./2.4 | field_loop.cpp:L329 |
| `vy0` | Real | 0.0 | turb_mhd_amr_wave.cpp:L57 |
| `vz` | Real | 0.0 | field_loop.cpp:L330 |
| `vz0` | Real | 0.0 | turb_mhd_amr_wave.cpp:L58 |
| `wave_amplitude` | Real | 0.0 | turb_mhd_amr_wave.cpp:L74 |
| `wave_direction` | int | 0 | turb_mhd_amr_wave.cpp:L75 |
| `wave_speed` | Real | 0.1 | turb_mhd_amr_wave.cpp:L72 |
| `wave_width` | Real | 0.2 | turb_mhd_amr_wave.cpp:L73 |
| `width` | Real | required | rad_snake.cpp:L193 |
| `x01` | Real | 3.0 | current_sheet.cpp:L44 |
| `x_CM` | Real | required | sgrid_bns.cpp:L154 |
| `x_center` | Real | 0.0 | sfb_turb.cpp:L61 |
| `xmax1` | Real | required | sgrid_bns.cpp:L157 |
| `xmax2` | Real | required | sgrid_bns.cpp:L158 |
| `y0` | Real | 0.0 | kh.cpp:L40 |
| `y1` | Real | 1.0 | kh.cpp:L41 |
| `y_center` | Real | 0.0 | sfb_turb.cpp:L62 |
| `z_center` | Real | 0.0 | sfb_turb.cpp:L63 |

## Input Block: `<radiation>`
**Used by**: rad_shadow.cpp, radiation.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `affect_fluid` | bool | true | radiation.cpp:L102 |
| `angular_fluxes` | bool | true | radiation.cpp:L115 |
| `arad` | Real | required | radiation.cpp:L100 |
| `beam_source` | bool | false | radiation.cpp:L109 |
| `compton` | bool | false | radiation.cpp:L89 |
| `fixed_fluid` | bool | false | radiation.cpp:L106 |
| `kappa_a` | Real | required | radiation.cpp:L86 |
| `kappa_p` | Real | required | radiation.cpp:L87 |
| `kappa_s` | Real | required | radiation.cpp:L83 |
| `n_0_floor` | Real | 0.1 | radiation.cpp:L116 |
| `nlevel` | int | required | rad_shadow.cpp:L61 |
| `power_opacity` | bool | false | radiation.cpp:L84 |
| `rad_source` | bool | true | radiation.cpp:L75 |
| `reconstruct` | string | plm | radiation.cpp:L167 |
| `rotate_geo` | bool | required | rad_shadow.cpp:L62 |

## Input Block: `<shearing_box>`
**Used by**: srcterms.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `omega0` | Real | required | srcterms.cpp:L103 |
| `qshear` | Real | required | srcterms.cpp:L102 |

## Input Block: `<time>`
**Used by**: build_tree.cpp, driver.cpp, hydro.cpp, part_static_turb.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `cfl_number` | Real | 0.8 | part_static_turb.cpp:L352 |
| `evolution` | string | required | hydro.cpp:L94 |
| `integrator` | string | rk2 | driver.cpp:L92 |
| `ndiag` | int | 1 | driver.cpp:L95 |
| `nlim` | int | -1 | driver.cpp:L94 |
| `start_time` | Real | 0.0 | build_tree.cpp:L304 |
| `tlim` | Real | required | driver.cpp:L93 |

## Input Block: `<turb_driving>`
**Used by**: turb_driver.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `basis_type` | int | 0 | turb_driver.cpp:L43 |
| `constant_edot` | bool | true | turb_driver.cpp:L103 |
| `dedt` | Real | 0.0 | turb_driver.cpp:L89 |
| `driving_type` | int | 0 | turb_driver.cpp:L76 |
| `dt_turb_thresh` | Real | 1e-6 | turb_driver.cpp:L95 |
| `dt_turb_update` | Real | 0.01 | turb_driver.cpp:L93 |
| `exp_prl` | Real | 0.0 | turb_driver.cpp:L87 |
| `exp_prp` | Real | 5.0/3.0 | turb_driver.cpp:L86 |
| `expo` | Real | 5.0/3.0 | turb_driver.cpp:L85 |
| `kpeak` | Real | 4.0*M_PI | turb_driver.cpp:L72 |
| `lmax` | int | 10 | turb_driver.cpp:L45 |
| `max_kx` | int | nhigh | turb_driver.cpp:L81 |
| `max_ky` | int | nhigh | turb_driver.cpp:L83 |
| `max_kz` | int | nhigh | turb_driver.cpp:L79 |
| `min_kx` | int | 0 | turb_driver.cpp:L80 |
| `min_ky` | int | 0 | turb_driver.cpp:L82 |
| `min_kz` | int | 0 | turb_driver.cpp:L78 |
| `nhigh` | int | 3 | turb_driver.cpp:L70 |
| `nlow` | int | 1 | turb_driver.cpp:L69 |
| `nmax` | int | 10 | turb_driver.cpp:L46 |
| `r0_turb` | Real | 0.5 | turb_driver.cpp:L47 |
| `rseed` | int | -1 | turb_driver.cpp:L100 |
| `sol_fraction` | Real | 1.0 | turb_driver.cpp:L97 |
| `spect_form` | int | 1 | turb_driver.cpp:L74 |
| `tcorr` | Real | 0.0 | turb_driver.cpp:L91 |
| `tdriv_duration` | Real | tcorr | turb_driver.cpp:L116 |
| `tdriv_start` | Real | 0. | turb_driver.cpp:L121 |
| `turb_flag` | int | 2 | turb_driver.cpp:L114 |
| `x_turb_center` | Real | 0.0 | turb_driver.cpp:L109 |
| `x_turb_scale_height` | Real | -1.0 | turb_driver.cpp:L106 |
| `y_turb_center` | Real | 0.0 | turb_driver.cpp:L110 |
| `y_turb_scale_height` | Real | -1.0 | turb_driver.cpp:L107 |
| `z_turb_center` | Real | 0.0 | turb_driver.cpp:L111 |
| `z_turb_scale_height` | Real | -1.0 | turb_driver.cpp:L108 |

## Input Block: `<units>`
**Used by**: units.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `bhmass_msun` | Real | required | units.cpp:L27 |
| `density_cgs` | Real | required | units.cpp:L26 |
| `length_cgs` | Real | 1.0 | units.cpp:L17 |
| `mass_cgs` | Real | 1.0 | units.cpp:L18 |
| `mu` | Real | 1.0 | units.cpp:L20 |
| `time_cgs` | Real | 1.0 | units.cpp:L19 |

## Input Block: `<z4c>`
**Used by**: compact_object_tracker.cpp, z4c.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `chi_div_floor` | Real | -1000.0 | z4c.cpp:L122 |
| `chi_psi_power` | Real | -4.0 | z4c.cpp:L121 |
| `damp_kappa1` | Real | 0.0 | z4c.cpp:L125 |
| `damp_kappa2` | Real | 0.0 | z4c.cpp:L126 |
| `diss` | Real | 0.0 | z4c.cpp:L123 |
| `eps_floor` | Real | 1e-12 | z4c.cpp:L124 |
| `extraction_nlev` | Real | 10 | z4c.cpp:L174 |
| `extrap_order` | int | 2 | z4c.cpp:L145 |
| `filename` | string | co_ | compact_object_tracker.cpp:L38 |
| `lapse_advect` | Real | 1.0 | z4c.cpp:L131 |
| `lapse_harmonic` | Real | 0.0 | z4c.cpp:L129 |
| `lapse_harmonicf` | Real | 1.0 | z4c.cpp:L128 |
| `lapse_oplog` | Real | 2.0 | z4c.cpp:L130 |
| `nrad_wave_extraction` | Real | 1 | z4c.cpp:L173 |
| `shift_Gamma` | Real | 1.0 | z4c.cpp:L132 |
| `shift_H` | Real | 0.0 | z4c.cpp:L136 |
| `shift_advect` | Real | 1.0 | z4c.cpp:L133 |
| `shift_alpha2Gamma` | Real | 0.0 | z4c.cpp:L135 |
| `shift_eta` | Real | 2.0 | z4c.cpp:L138 |
| `use_z4c` | bool | true | z4c.cpp:L140 |
| `user_Sbc` | bool | false | z4c.cpp:L142 |
| `waveform_dt` | Real | 1 | z4c.cpp:L182 |

## Input Block: `<z4c_amr>`
**Used by**: z4c_amr.cpp

| Parameter | Type | Default | Source |
|-----------|------|---------|--------|
| `chi_min` | Real | 0.2 | z4c_amr.cpp:L35 |
| `dchi_max` | Real | 0.1 | z4c_amr.cpp:L38 |
| `method` | string | trivial | z4c_amr.cpp:L28 |


## Summary
- Total input blocks: 19
- Total parameters: 340
