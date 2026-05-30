//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.cpp
//! \brief implementation of Particles class constructor and assorted other functions

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "particle_restart.hpp"
#include "particles.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace particles {
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Particles::Particles(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack),
    dtnew(std::numeric_limits<float>::max()),
    relativistic_restart_v2_loaded(false),
    relativistic_restart_mesh_cycle(0),
    relativistic_restart_mesh_time(0.0),
    relativistic_restart_mesh_dt(0.0),
    relativistic_restart_config_fingerprint(0),
    legacy_restart_mesh_time(0.0),
    legacy_restart_mesh_dt(0.0),
    check_consistency(false),
    check_motion_bounds(false),
    log_performance(false),
    validate_amr_lookup(false),
    subcycle(false),
    subcycle_strict(true),
    amr_lookup_max_cells(20000000),
    subcycle_max_steps(8),
    subcycle_cell_fraction(0.5),
    subcycle_meshblock_fraction(0.5),
    subcycle_gyro_fraction(0.25),
    subcycle_electric_kick_max(0.1),
    c_model(0.0),
    alpha_s(0.0),
    relativistic_timestep_bound_dirty(true),
    relativistic_timestep_bound_override(false),
    last_subcycle_steps(1),
    last_push_dt(0.0),
    consistency_mode(ParticlesConsistencyMode::none),
    amr_remap_mode(ParticlesAMRRemapMode::device_table),
    exchange_mode(ParticlesExchangeMode::alltoall_counts),
    interpolation(ParticleInterpolation::tsc),
    relativistic_field_source(RelativisticFieldSource::prescribed_test),
    relativistic_temporal_sampling(RelativisticTemporalSampling::none),
    reference_counts_set(false),
    reference_nprtcl_total(0) {
  // check this is at least a 2D problem
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle module only works in 2D/3D" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // read number of particles per cell, and calculate number of particles this pack
  Real ppc = pin->GetOrAddReal("particles","ppc",1.0);
  nspecies = pin->GetOrAddInteger("particles","nspecies",1);
  if (nspecies < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<particles>/nspecies must be >= 1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string evolution_t = pin->GetOrAddString("time","evolution","dynamic");
  is_dynamic = (evolution_t.compare("static") != 0);

  // compute number of particles as real number, since ppc can be < 1
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
  Real r_npart = ppc*static_cast<Real>((pmy_pack->nmb_thispack)*ncells);
  // then cast to integer
  nprtcl_perspec_thispack = static_cast<int>(r_npart);
  nprtcl_thispack = nprtcl_perspec_thispack*nspecies;

  prtcl_rst_flag = pin->GetOrAddInteger("problem","prtcl_rst_flag",0);
  bool legacy_check_consistency = pin->GetOrAddBoolean("particles",
                                                       "check_consistency",false);
  std::string check_mode = pin->GetOrAddString("particles","check_consistency_mode",
                                               legacy_check_consistency ?
                                               "full" : "none");
  if (check_mode.compare("none") == 0 || check_mode.compare("false") == 0) {
    consistency_mode = ParticlesConsistencyMode::none;
  } else if (check_mode.compare("counts") == 0) {
    consistency_mode = ParticlesConsistencyMode::counts;
  } else if (check_mode.compare("local") == 0) {
    consistency_mode = ParticlesConsistencyMode::local;
  } else if (check_mode.compare("full") == 0 || check_mode.compare("true") == 0) {
    consistency_mode = ParticlesConsistencyMode::full;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<particles>/check_consistency_mode = '"
              << check_mode << "' not recognized" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  check_consistency = (consistency_mode != ParticlesConsistencyMode::none);
  check_motion_bounds = pin->GetOrAddBoolean("particles","check_motion_bounds",false);
  log_performance = pin->GetOrAddBoolean("particles","log_performance",false);
  validate_amr_lookup = pin->GetOrAddBoolean("particles","validate_amr_lookup",false);
  amr_lookup_max_cells = pin->GetOrAddInteger("particles","amr_lookup_max_cells",
                                              20000000);
  subcycle = pin->GetOrAddBoolean("particles","subcycle",false);
  subcycle_max_steps = pin->GetOrAddInteger("particles","subcycle_max_steps",8);
  subcycle_cell_fraction = pin->GetOrAddReal("particles","subcycle_cell_fraction",0.5);
  subcycle_meshblock_fraction =
      pin->GetOrAddReal("particles","subcycle_meshblock_fraction",0.5);
  subcycle_gyro_fraction = pin->GetOrAddReal("particles","subcycle_gyro_fraction",0.25);
  bool has_relativistic_electric_kick_max =
      pin->DoesParameterExist("particles","subcycle_electric_kick_max");
  bool has_deferred_momentum_fraction =
      pin->DoesParameterExist("particles","subcycle_momentum_fraction");
  bool has_deferred_kinetic_fraction =
      pin->DoesParameterExist("particles","subcycle_kinetic_energy_fraction");
  subcycle_electric_kick_max =
      pin->GetOrAddReal("particles","subcycle_electric_kick_max",0.1);
  subcycle_strict = pin->GetOrAddBoolean("particles","subcycle_strict",true);
  if (subcycle_max_steps < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<particles>/subcycle_max_steps must be >= 1"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (subcycle_cell_fraction <= 0.0 || subcycle_meshblock_fraction <= 0.0 ||
      subcycle_gyro_fraction <= 0.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Particle subcycle fractions must be positive"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string amr_remap = pin->GetOrAddString("particles","amr_remap","device_table");
  if (amr_remap.compare("device_table") == 0) {
    amr_remap_mode = ParticlesAMRRemapMode::device_table;
  } else if (amr_remap.compare("host_tree") == 0) {
    amr_remap_mode = ParticlesAMRRemapMode::host_tree;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<particles>/amr_remap = '" << amr_remap
              << "' not recognized" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string exchange = pin->GetOrAddString("particles","exchange_mode",
                                             "alltoall_counts");
  if (exchange.compare("alltoall_counts") == 0) {
    exchange_mode = ParticlesExchangeMode::alltoall_counts;
  } else if (exchange.compare("allgather") == 0) {
    exchange_mode = ParticlesExchangeMode::allgather;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<particles>/exchange_mode = '" << exchange
              << "' not recognized" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string ppush = pin->GetString("particles","pusher");
  bool has_relativistic_key =
      pin->DoesParameterExist("particles","relativistic_field_source") ||
      pin->DoesParameterExist("particles","relativistic_temporal_sampling") ||
      pin->DoesParameterExist("particles","relativistic_mhd_temporal_mode") ||
      pin->DoesParameterExist("particles","c_model") ||
      pin->DoesParameterExist("particles","alpha_s") ||
      has_relativistic_electric_kick_max ||
      has_deferred_momentum_fraction ||
      has_deferred_kinetic_fraction ||
      pin->DoesParameterExist("particles","relativistic_background");
  bool has_relativistic_problem_key =
      pin->DoesParameterExist("problem","relativistic_initial_state") ||
      pin->DoesParameterExist("problem","w0x") ||
      pin->DoesParameterExist("problem","w0y") ||
      pin->DoesParameterExist("problem","w0z") ||
      pin->DoesParameterExist("problem","cE0x") ||
      pin->DoesParameterExist("problem","cE0y") ||
      pin->DoesParameterExist("problem","cE0z") ||
      pin->DoesParameterExist("problem","relativistic_gather_diagnostic_only");
  if (ppush.compare("relativistic_hc") != 0 &&
      (has_relativistic_key || has_relativistic_problem_key)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Relativistic-only input keys require "
              << "pusher = relativistic_hc" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ppush.compare("relativistic_hc") == 0) {
    int qualified_prst = 0;
    int qualified_rst = 0;
    for (auto it = pin->block.begin(); it != pin->block.end(); ++it) {
      if (it->block_name.compare(0, 6, "output") == 0) {
        std::string file_type = pin->GetString(it->block_name, "file_type");
        if (file_type.compare("prst") == 0) {
          ++qualified_prst;
          continue;
        }
        if (file_type.compare("rst") == 0) {
          ++qualified_rst;
          continue;
        }
        if (file_type.compare("ppd") == 0 ||
            file_type.compare("df") == 0 ||
            file_type.compare("dxh") == 0 ||
            file_type.compare("drh") == 0 ||
            file_type.compare("dparh") == 0 ||
            file_type.compare("pmom") == 0 ||
            file_type.compare("pspec") == 0 ||
            file_type.compare("pspec2") == 0 ||
            file_type.compare("psamp") == 0) {
          continue;
        }
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' staged "
                  << "diagnostic contract rejects <" << it->block_name
                  << "> file_type = " << file_type
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (qualified_prst != qualified_rst || qualified_prst > 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher 'relativistic_hc' requires exactly "
                << "one paired <output>/file_type = prst and rst block, or neither"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (prtcl_rst_flag) {
    std::string prst_fname = pin->GetString("problem","prtcl_res_file");
    if (ppush.compare("relativistic_hc") == 0) {
      try {
        std::string mesh_restart =
            pin->GetString("problem","relativistic_paired_mesh_restart_file");
        restart::V2Shard shard = restart::ReadV2Shard(
            prst_fname, mesh_restart, global_variable::my_rank,
            global_variable::nranks,
            restart::ComputeTopologyHash(pmy_pack->pmesh->rank_eachmb,
                                         pmy_pack->pmesh->nmb_total));
        nprtcl_thispack = static_cast<int>(shard.header.local_count);
        nprtcl_perspec_thispack = nprtcl_thispack/nspecies;
        relativistic_restart_v2_loaded = true;
        relativistic_restart_mesh_cycle = shard.header.mesh_cycle;
        relativistic_restart_mesh_time = shard.header.mesh_time;
        relativistic_restart_mesh_dt = shard.header.mesh_dt;
        relativistic_restart_config_fingerprint = shard.header.config_fingerprint;
        dtnew = shard.header.particle_dtnew;
        relativistic_restart_idata = std::move(shard.idata);
        relativistic_restart_rdata = std::move(shard.rdata);
      } catch (const std::exception &error) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Relativistic typed-v2 particle restart rejected: "
                  << error.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      try {
        restart::LegacyV1Shard<Real> shard =
            restart::ReadLegacyV1Shard<Real>(prst_fname, global_variable::my_rank);
        nprtcl_thispack = static_cast<int>(shard.rdata.size()/17);
        nprtcl_perspec_thispack = nprtcl_thispack/nspecies;
        legacy_restart_mesh_time = shard.mesh_time;
        legacy_restart_mesh_dt = shard.mesh_dt;
        legacy_restart_rdata = std::move(shard.rdata);
      } catch (const std::exception &error) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Legacy particle restart rejected: "
                  << error.what() << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }

  // select particle type
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("cosmic_ray") == 0) {
      particle_type = ParticleType::cosmic_ray;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle type = '" << ptype << "' not recognized"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // select pusher algorithm
  {
    std::string interp = pin->GetOrAddString("particles","interpolation", "tsc");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
    } else if (ppush.compare("boris") == 0) {
      if (pmy_pack->pmhd == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'boris' requires an <mhd> block"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (interp.compare("tsc") == 0) {
        interpolation = ParticleInterpolation::tsc;
      } else if (interp.compare("lin") == 0 || interp.compare("lin_legacy") == 0) {
        interpolation = ParticleInterpolation::lin_legacy;
      } else if (interp.compare("trilinear") == 0) {
        interpolation = ParticleInterpolation::trilinear;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Interpolation method = '" << interp
                  << "' not recognized" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      pusher = ParticlesPusher::boris;
    } else if (ppush.compare("relativistic_hc") == 0) {
      bool gather_diagnostic_only =
          pin->GetOrAddBoolean("problem","relativistic_gather_diagnostic_only",false);
      if (gather_diagnostic_only &&
          (pin->GetString("problem","pgen_name").compare(
               "cr_relativistic_runtime_gather_test") != 0 ||
           pin->GetInteger("time","nlim") != 0)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<problem>/relativistic_gather_diagnostic_only "
                  << "requires pgen_name = cr_relativistic_runtime_gather_test and "
                  << "<time>/nlim = 0" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pin->DoesParameterExist("particles","relativistic_background")) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<particles>/relativistic_background is unsupported; "
                  << "select <particles>/relativistic_field_source explicitly"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pin->DoesParameterExist("particles","relativistic_mhd_temporal_mode")) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<particles>/relativistic_mhd_temporal_mode is "
                  << "unsupported; select <particles>/relativistic_temporal_sampling"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pmy_pack->pmhd == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' requires an "
                  << "<mhd> block" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (nspecies != 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' currently requires "
                  << "<particles>/nspecies = 1" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!gather_diagnostic_only &&
          (global_variable::nranks != 1 || pmy_pack->pmesh->multilevel)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' execution currently "
                  << "requires a serial uniform-level mesh; "
                  << "MPI, SMR, and AMR remain closed until qualification" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!(pmy_pack->pmesh->three_d) || !(pmy_pack->pmesh->strictly_periodic)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' currently requires "
                  << "a 3D strictly periodic mesh" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pmy_pack->pcoord->is_special_relativistic ||
          pmy_pack->pcoord->is_general_relativistic ||
          pmy_pack->pcoord->is_dynamical_relativistic) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' currently supports "
                  << "Newtonian coordinates only" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      bool coupled_manufactured_test =
          pin->GetOrAddBoolean("problem","relativistic_coupled_manufactured_test",false);
      bool allow_manufactured_user_src =
          coupled_manufactured_test &&
          std::string(PROBLEM_GENERATOR).compare(
              "unit_tests/cr_relativistic_coupled_runtime_test") == 0 &&
          pin->GetString("problem","pgen_name").compare(
              "cr_relativistic_coupled_runtime_test") == 0;
      if (coupled_manufactured_test && !allow_manufactured_user_src) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<problem>/relativistic_coupled_manufactured_test "
                  << "requires the compiled coupled-runtime test harness and "
                  << "pgen_name = cr_relativistic_coupled_runtime_test"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("radiation") ||
          pin->DoesBlockExist("ion-neutral") || pin->DoesBlockExist("turb_driving") ||
          pin->DoesBlockExist("mhd_srcterms") || pin->DoesBlockExist("hydro_srcterms") ||
          pin->DoesBlockExist("rad_srcterms") || pin->DoesBlockExist("shearing_box") ||
          pin->DoesBlockExist("adm") || pin->DoesBlockExist("z4c") ||
          (pin->GetOrAddBoolean("problem","user_srcs",false) &&
           !allow_manufactured_user_src)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' contract rejects "
                  << "unqualified coupled physics blocks" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pin->DoesParameterExist("mhd","viscosity") ||
          pin->DoesParameterExist("mhd","ohmic_resistivity") ||
          pin->DoesParameterExist("mhd","conductivity") ||
          pin->DoesParameterExist("mhd","tdep_conductivity")) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' contract rejects "
                  << "unqualified MHD diffusion" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (pin->DoesParameterExist("mhd","nscalars") &&
          pin->GetInteger("mhd","nscalars") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' contract rejects "
                  << "unqualified MHD passive scalars" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::string eos = pin->GetString("mhd","eos");
      if (eos.compare("ideal") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' requires "
                  << "<mhd>/eos = ideal" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (subcycle && !subcycle_strict) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' requires "
                  << "<particles>/subcycle_strict = true; cap clipping is unsupported"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (has_deferred_momentum_fraction || has_deferred_kinetic_fraction) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' does not support "
                  << "<particles>/subcycle_momentum_fraction or "
                  << "<particles>/subcycle_kinetic_energy_fraction; both require a "
                  << "separately reviewed nonsingular near-rest scale" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!(pin->DoesParameterExist("particles","interpolation")) ||
          interp.compare("trilinear") != 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' requires explicit "
                  << "<particles>/interpolation = trilinear" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::string field_source =
          pin->GetOrAddString("particles","relativistic_field_source","");
      if (field_source.compare("prescribed_test") == 0) {
        relativistic_field_source = RelativisticFieldSource::prescribed_test;
        if (coupled_manufactured_test) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<problem>/relativistic_coupled_manufactured_test "
                    << "requires relativistic_field_source = mhd_ideal" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (pin->DoesParameterExist("particles","relativistic_temporal_sampling")) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<particles>/relativistic_temporal_sampling applies "
                    << "only to relativistic_field_source = mhd_ideal" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (evolution_t.compare("kinematic") != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "Particle pusher 'relativistic_hc' prescribed-test "
                    << "contract requires <time>/evolution = kinematic" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else if (field_source.compare("mhd_ideal") == 0) {
        relativistic_field_source = RelativisticFieldSource::mhd_ideal;
        if (pmy_pack->nmb_thispack != 1) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "relativistic_field_source = mhd_ideal currently "
                    << "requires exactly one MeshBlock until same-level migration is "
                    << "qualified in Phase 8" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (gather_diagnostic_only) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<problem>/relativistic_gather_diagnostic_only is a "
                    << "prescribed-test-only Phase-3 exception" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (pin->DoesParameterExist("problem","cE0x") ||
            pin->DoesParameterExist("problem","cE0y") ||
            pin->DoesParameterExist("problem","cE0z")) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "relativistic_field_source = mhd_ideal rejects "
                    << "prescribed-test-only cE0x, cE0y, and cE0z" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (!(pin->DoesParameterExist("particles","relativistic_temporal_sampling"))) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "relativistic_field_source = mhd_ideal requires "
                    << "<particles>/relativistic_temporal_sampling = frozen_tn"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        std::string temporal_mode =
            pin->GetString("particles","relativistic_temporal_sampling");
        if (temporal_mode.compare("frozen_tn") != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<particles>/relativistic_temporal_sampling = '"
                    << temporal_mode << "' not recognized; expected frozen_tn"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        relativistic_temporal_sampling = RelativisticTemporalSampling::frozen_tn;
        if (evolution_t.compare("dynamic") != 0) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "relativistic_field_source = mhd_ideal requires "
                    << "<time>/evolution = dynamic" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<particles>/relativistic_field_source = '"
                  << field_source << "' not recognized; expected prescribed_test or "
                  << "mhd_ideal"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!(pin->DoesParameterExist("particles","c_model"))) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' requires "
                  << "<particles>/c_model" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      c_model = pin->GetReal("particles","c_model");
      if (!std::isfinite(c_model) || c_model <= 0.0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<particles>/c_model must be finite and > 0 for "
                  << "particle pusher 'relativistic_hc'" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (relativistic_field_source == RelativisticFieldSource::mhd_ideal &&
          c_model != 1.0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "relativistic_field_source = mhd_ideal currently "
                  << "requires <particles>/c_model = 1.0 until non-unit coupled "
                  << "normalization is separately qualified" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!(pin->DoesParameterExist("particles","alpha_s"))) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' requires "
                  << "<particles>/alpha_s" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      alpha_s = pin->GetReal("particles","alpha_s");
      if (!std::isfinite(alpha_s) || alpha_s == 0.0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<particles>/alpha_s must be finite and nonzero for "
                  << "particle pusher 'relativistic_hc'" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::uint32_t field_source_id =
          (relativistic_field_source == RelativisticFieldSource::prescribed_test) ?
          1U : 2U;
      std::uint32_t temporal_sampling_id =
          (relativistic_temporal_sampling == RelativisticTemporalSampling::none) ?
          0U : 1U;
      if (relativistic_restart_v2_loaded &&
          relativistic_restart_config_fingerprint !=
          restart::ComputeRelativisticConfigFingerprint(
              c_model, alpha_s, field_source_id, temporal_sampling_id,
              subcycle, subcycle_strict, subcycle_max_steps,
              subcycle_cell_fraction, subcycle_meshblock_fraction,
              subcycle_gyro_fraction, subcycle_electric_kick_max)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Relativistic typed-v2 particle restart semantic "
                  << "configuration does not match the active input" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!std::isfinite(subcycle_cell_fraction) ||
          !std::isfinite(subcycle_meshblock_fraction) ||
          !std::isfinite(subcycle_gyro_fraction) ||
          !std::isfinite(subcycle_electric_kick_max) ||
          subcycle_cell_fraction > 0.5 || subcycle_meshblock_fraction > 0.5 ||
          subcycle_electric_kick_max <= 0.0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Particle pusher 'relativistic_hc' requires finite "
                  << "subcycle fractions, cell and MeshBlock fractions <= 0.5, and "
                  << "<particles>/subcycle_electric_kick_max > 0" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      pusher = ParticlesPusher::relativistic_hc;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher must be specified in <particles> block"
                <<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // set dimensions of particle arrays. Note particles only work in 2D/3D
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particles only work in 2D/3D, but 1D problem initialized" <<std::endl;
    std::exit(EXIT_FAILURE);
  }
  switch (particle_type) {
    case ParticleType::cosmic_ray:
      {
        nrdata = (pusher == ParticlesPusher::relativistic_hc) ?
            static_cast<int>(restart::kV2RealFields) : 14;
        nidata = (pusher == ParticlesPusher::relativistic_hc) ?
            static_cast<int>(restart::kV2IntegerFields) : 3;
        break;
      }
    default:
      break;
  }
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  // allocate boundary object
  pbval_part = new ParticlesBoundaryValues(this, pin);
}

//----------------------------------------------------------------------------------------
// destructor

Particles::~Particles() {
}

//----------------------------------------------------------------------------------------
// CreateParticleTags()
// Assigns tags to particles (unique integer).  Note that tracked particles are always
// those with tag numbers less than ntrack.

void Particles::CreateParticleTags(ParameterInput *pin) {
  if (prtcl_rst_flag) {return;}

  std::string assign = pin->GetOrAddString("particles","assign_tag","index_order");
  int nprtcl_total_perspec = pmy_pack->pmesh->nprtcl_total/nspecies;

  // tags are assigned sequentially within this rank, starting at 0 with rank=0
  if (assign.compare("index_order") == 0) {
    int tagstart = 0;
    for (int n=1; n<=global_variable::my_rank; ++n) {
      tagstart += pmy_pack->pmesh->nprtcl_eachrank[n-1]/nspecies;
    }

    auto &pi = prtcl_idata;
    int nper = nprtcl_perspec_thispack;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      int spec = (nper > 0) ? p/nper : 0;
      int p_in_spec = (nper > 0) ? p - spec*nper : p;
      pi(PTAG,p) = spec*nprtcl_total_perspec + tagstart + p_in_spec;
    });

  // tags are assigned sequentially across ranks
  } else if (assign.compare("rank_order") == 0) {
    int myrank = global_variable::my_rank;
    int nranks = global_variable::nranks;
    auto &pi = prtcl_idata;
    int nper = nprtcl_perspec_thispack;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      int spec = (nper > 0) ? p/nper : 0;
      int p_in_spec = (nper > 0) ? p - spec*nper : p;
      pi(PTAG,p) = spec*nprtcl_total_perspec + myrank + nranks*p_in_spec;
    });

  // tag algorithm not recognized, so quit with error
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle tag assignment type = '" << assign << "' not recognized"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

namespace {

void WrapCoordinate(Real &x, Real xmin, Real xmax) {
  Real len = xmax - xmin;
  if (len <= 0.0) {return;}
  while (x < xmin) {x += len;}
  while (x >= xmax) {x -= len;}
}

KOKKOS_INLINE_FUNCTION
void WrapCoordinateDevice(Real &x, Real xmin, Real xmax) {
  Real len = xmax - xmin;
  if (len <= 0.0) {return;}
  while (x < xmin) {x += len;}
  while (x >= xmax) {x -= len;}
}

KOKKOS_INLINE_FUNCTION
int LocationIndexDevice(Real x, Real xmin, Real xmax, int nloc) {
  Real frac = (x - xmin)/(xmax - xmin);
  std::int64_t ix = static_cast<std::int64_t>(frac*static_cast<Real>(nloc));
  ix = (ix < 0) ? 0 : ix;
  ix = (ix > nloc - 1) ? nloc - 1 : ix;
  return static_cast<int>(ix);
}

void LogicalLocationToRegion(Mesh *pm, const LogicalLocation &lloc, RegionSize &rs) {
  RegionSize &ms = pm->mesh_size;
  int lev_offset = lloc.level - pm->root_level;

  std::int32_t nmbx1 = pm->nmb_rootx1 << lev_offset;
  rs.x1min = (lloc.lx1 == 0) ? ms.x1min : LeftEdgeX(lloc.lx1, nmbx1,
                                                     ms.x1min, ms.x1max);
  rs.x1max = (lloc.lx1 == nmbx1 - 1) ? ms.x1max : LeftEdgeX(lloc.lx1 + 1, nmbx1,
                                                            ms.x1min, ms.x1max);

  if (pm->multi_d) {
    std::int32_t nmbx2 = pm->nmb_rootx2 << lev_offset;
    rs.x2min = (lloc.lx2 == 0) ? ms.x2min : LeftEdgeX(lloc.lx2, nmbx2,
                                                       ms.x2min, ms.x2max);
    rs.x2max = (lloc.lx2 == nmbx2 - 1) ? ms.x2max : LeftEdgeX(lloc.lx2 + 1, nmbx2,
                                                              ms.x2min, ms.x2max);
  } else {
    rs.x2min = ms.x2min;
    rs.x2max = ms.x2max;
  }

  if (pm->three_d) {
    std::int32_t nmbx3 = pm->nmb_rootx3 << lev_offset;
    rs.x3min = (lloc.lx3 == 0) ? ms.x3min : LeftEdgeX(lloc.lx3, nmbx3,
                                                       ms.x3min, ms.x3max);
    rs.x3max = (lloc.lx3 == nmbx3 - 1) ? ms.x3max : LeftEdgeX(lloc.lx3 + 1, nmbx3,
                                                              ms.x3min, ms.x3max);
  } else {
    rs.x3min = ms.x3min;
    rs.x3max = ms.x3max;
  }
}

bool ContainsParticle(Mesh *pm, const RegionSize &rs, const LogicalLocation &lloc,
                      Real x1, Real x2, Real x3) {
  int lev_offset = lloc.level - pm->root_level;
  std::int32_t nmbx1 = pm->nmb_rootx1 << lev_offset;
  std::int32_t nmbx2 = pm->nmb_rootx2 << lev_offset;
  std::int32_t nmbx3 = pm->nmb_rootx3 << lev_offset;

  bool in_x1 = (x1 >= rs.x1min) && (x1 < rs.x1max ||
               (lloc.lx1 == nmbx1 - 1 && x1 <= rs.x1max));
  bool in_x2 = (!pm->multi_d) || ((x2 >= rs.x2min) && (x2 < rs.x2max ||
               (lloc.lx2 == nmbx2 - 1 && x2 <= rs.x2max)));
  bool in_x3 = (!pm->three_d) || ((x3 >= rs.x3min) && (x3 < rs.x3max ||
               (lloc.lx3 == nmbx3 - 1 && x3 <= rs.x3max)));
  return in_x1 && in_x2 && in_x3;
}

void FatalConsistencyError(const std::string &label, const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Particle consistency check '" << label
            << "' failed: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

void FatalMotionError(const std::string &label, const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Particle motion-bounds check '" << label
            << "' failed: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

void FatalAMRRemapError(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << "Particle AMR remap failed: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

bool BuildAMRLookupTable(Mesh *pm, int max_entries, DvceArray1D<int> &gid_table,
                         DvceArray1D<int> &rank_table, int &nloc1, int &nloc2,
                         int &nloc3) {
  int lev_offset = pm->max_level - pm->root_level;
  nloc1 = pm->nmb_rootx1 << lev_offset;
  nloc2 = pm->multi_d ? (pm->nmb_rootx2 << lev_offset) : 1;
  nloc3 = pm->three_d ? (pm->nmb_rootx3 << lev_offset) : 1;

  int64_t ntable = static_cast<int64_t>(nloc1)*static_cast<int64_t>(nloc2)*
                   static_cast<int64_t>(nloc3);
  if (ntable > static_cast<int64_t>(max_entries)) {return false;}

  HostArray1D<int> gid_h("amr_lookup_gid_h", ntable);
  HostArray1D<int> rank_h("amr_lookup_rank_h", ntable);
  Kokkos::deep_copy(gid_h, -1);
  Kokkos::deep_copy(rank_h, -1);

  for (int gid=0; gid<pm->nmb_total; ++gid) {
    LogicalLocation &lloc = pm->lloc_eachmb[gid];
    int shift = pm->max_level - lloc.level;
    int ilo = lloc.lx1 << shift;
    int ihi = (lloc.lx1 + 1) << shift;
    int jlo = pm->multi_d ? (lloc.lx2 << shift) : 0;
    int jhi = pm->multi_d ? ((lloc.lx2 + 1) << shift) : 1;
    int klo = pm->three_d ? (lloc.lx3 << shift) : 0;
    int khi = pm->three_d ? ((lloc.lx3 + 1) << shift) : 1;
    for (int k=klo; k<khi; ++k) {
      for (int j=jlo; j<jhi; ++j) {
        for (int i=ilo; i<ihi; ++i) {
          int64_t idx = static_cast<int64_t>(i) +
                        static_cast<int64_t>(nloc1)*
                        (static_cast<int64_t>(j) +
                         static_cast<int64_t>(nloc2)*static_cast<int64_t>(k));
          gid_h(idx) = gid;
          rank_h(idx) = pm->rank_eachmb[gid];
        }
      }
    }
  }

  for (int64_t n=0; n<ntable; ++n) {
    if (gid_h(n) < 0 || rank_h(n) < 0) {
      FatalAMRRemapError("device lookup table has an unfilled logical slot");
    }
  }

  gid_table = DvceArray1D<int>("amr_lookup_gid", ntable);
  rank_table = DvceArray1D<int>("amr_lookup_rank", ntable);
  Kokkos::deep_copy(gid_table, gid_h);
  Kokkos::deep_copy(rank_table, rank_h);
  return true;
}

} // namespace

//----------------------------------------------------------------------------------------
// SetConsistencyReference()
// Store the expected total and per-species particle counts for debug checks.

void Particles::SetConsistencyReference() {
  Mesh *pm = pmy_pack->pmesh;
  pm->UpdateParticleCounts();
  reference_nprtcl_eachspec.assign(nspecies, 0);
  auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  std::vector<int> local_counts(nspecies, 0);
  for (int p=0; p<nprtcl_thispack; ++p) {
    int spec = pi_h(PSP,p);
    if (spec >= 0 && spec < nspecies) {local_counts[spec]++;}
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(local_counts.data(), reference_nprtcl_eachspec.data(), nspecies,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  reference_nprtcl_eachspec = local_counts;
#endif
  reference_nprtcl_total = pm->nprtcl_total;
  reference_counts_set = true;
}

//----------------------------------------------------------------------------------------
// CheckConsistency()
// Host-side debug checks for particle ownership, positions, tags, and conserved counts.

void Particles::CheckConsistency(const std::string &label, bool check_rank) {
  if (!check_consistency) {return;}

  Mesh *pm = pmy_pack->pmesh;
  pm->UpdateParticleCounts();
  auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  bool check_local_state = (consistency_mode == ParticlesConsistencyMode::local ||
                            consistency_mode == ParticlesConsistencyMode::full);
  HostArray2D<Real> pr_h;
  if (check_local_state) {
    pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  }
  std::vector<int> local_counts(nspecies, 0);
  std::vector<int64_t> local_keys;
  if (check_local_state) {local_keys.reserve(std::max(0, nprtcl_thispack));}

  int invalid_gid_errors = 0;
  int wrong_rank_errors = 0;
  int invalid_species_errors = 0;
  int nonfinite_errors = 0;
  int bounds_errors = 0;
  int local_duplicate_errors = 0;
  int count_errors = 0;
  int global_duplicate_errors = 0;
  std::string first_local_error;
  RegionSize rs;
  int myrank = global_variable::my_rank;
  for (int p=0; p<nprtcl_thispack; ++p) {
    int gid = pi_h(PGID,p);
    int tag = pi_h(PTAG,p);
    int spec = pi_h(PSP,p);
    if (gid < 0 || gid >= pm->nmb_total) {
      invalid_gid_errors++;
      if (first_local_error.empty()) {
        first_local_error = "invalid_gid rank=" + std::to_string(myrank) +
                            " p=" + std::to_string(p) +
                            " gid=" + std::to_string(gid);
      }
      continue;
    }
    if (check_rank && pm->rank_eachmb[gid] != myrank) {
      wrong_rank_errors++;
      if (first_local_error.empty()) {
        first_local_error = "wrong_rank rank=" + std::to_string(myrank) +
                            " p=" + std::to_string(p) +
                            " gid=" + std::to_string(gid) +
                            " owner=" + std::to_string(pm->rank_eachmb[gid]);
      }
    }
    if (spec < 0 || spec >= nspecies) {
      invalid_species_errors++;
      if (first_local_error.empty()) {
        first_local_error = "invalid_species rank=" + std::to_string(myrank) +
                            " p=" + std::to_string(p) +
                            " species=" + std::to_string(spec);
      }
    } else {
      local_counts[spec]++;
      if (check_local_state) {
        int64_t key = (static_cast<int64_t>(spec) << 32) ^
                      static_cast<unsigned int>(tag);
        local_keys.push_back(key);
      }
    }
    if (check_local_state) {
      for (int n=0; n<nrdata; ++n) {
        if (!std::isfinite(pr_h(n,p))) {
          nonfinite_errors++;
          if (first_local_error.empty()) {
            first_local_error = "nonfinite rank=" + std::to_string(myrank) +
                                " p=" + std::to_string(p) +
                                " field=" + std::to_string(n);
          }
        }
      }
      LogicalLocationToRegion(pm, pm->lloc_eachmb[gid], rs);
      if (!ContainsParticle(pm, rs, pm->lloc_eachmb[gid],
                            pr_h(IPX,p), pr_h(IPY,p), pr_h(IPZ,p))) {
        bounds_errors++;
        if (first_local_error.empty()) {
          first_local_error = "out_of_bounds rank=" + std::to_string(myrank) +
                              " p=" + std::to_string(p) +
                              " gid=" + std::to_string(gid) +
                              " x=(" + std::to_string(pr_h(IPX,p)) + "," +
                              std::to_string(pr_h(IPY,p)) + "," +
                              std::to_string(pr_h(IPZ,p)) + ") bounds=(" +
                              std::to_string(rs.x1min) + "," +
                              std::to_string(rs.x1max) + "; " +
                              std::to_string(rs.x2min) + "," +
                              std::to_string(rs.x2max) + "; " +
                              std::to_string(rs.x3min) + "," +
                              std::to_string(rs.x3max) + ")";
        }
      }
    }
  }

  std::vector<int> global_counts(nspecies, 0);
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(local_counts.data(), global_counts.data(), nspecies,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  global_counts = local_counts;
#endif
  if (reference_counts_set) {
    if (pm->nprtcl_total != reference_nprtcl_total) {count_errors++;}
    for (int spec=0; spec<nspecies; ++spec) {
      if (global_counts[spec] != reference_nprtcl_eachspec[spec]) {count_errors++;}
    }
  }

  if (check_local_state) {
    std::sort(local_keys.begin(), local_keys.end());
    for (std::size_t n=1; n<local_keys.size(); ++n) {
      if (local_keys[n] == local_keys[n-1]) {
        local_duplicate_errors++;
        if (first_local_error.empty()) {
          first_local_error = "local_duplicate_tag rank=" + std::to_string(myrank);
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  if (consistency_mode == ParticlesConsistencyMode::full) {
    int local_nkeys = static_cast<int>(local_keys.size());
    std::vector<int> recv_counts(global_variable::nranks, 0);
    MPI_Gather(&local_nkeys, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0,
               MPI_COMM_WORLD);
    std::vector<int> displs(global_variable::nranks, 0);
    std::vector<int64_t> all_keys;
    if (global_variable::my_rank == 0) {
      for (int n=1; n<global_variable::nranks; ++n) {
        displs[n] = displs[n-1] + recv_counts[n-1];
      }
      all_keys.resize(displs.back() + recv_counts.back());
    }
    MPI_Gatherv(local_keys.data(), local_nkeys, MPI_INT64_T,
                all_keys.data(), recv_counts.data(), displs.data(), MPI_INT64_T,
                0, MPI_COMM_WORLD);
    int duplicate_errors = 0;
    if (global_variable::my_rank == 0) {
      std::sort(all_keys.begin(), all_keys.end());
      for (std::size_t n=1; n<all_keys.size(); ++n) {
        if (all_keys[n] == all_keys[n-1]) {duplicate_errors++;}
      }
    }
    MPI_Bcast(&duplicate_errors, 1, MPI_INT, 0, MPI_COMM_WORLD);
    global_duplicate_errors += duplicate_errors;
  }
#endif

  int errors = invalid_gid_errors + wrong_rank_errors + invalid_species_errors +
               nonfinite_errors + bounds_errors + local_duplicate_errors +
               count_errors + global_duplicate_errors;
#if MPI_PARALLEL_ENABLED
  int global_errors = 0;
  MPI_Allreduce(&errors, &global_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &invalid_gid_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &wrong_rank_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &invalid_species_errors, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &nonfinite_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &bounds_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &local_duplicate_errors, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &count_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &global_duplicate_errors, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  errors = global_errors;
#endif
  if (errors > 0) {
    std::string details = std::to_string(errors) + " invariant violations";
    auto add_detail = [&](const std::string &name, int count) {
      if (count > 0) {
        details += "; " + name + "=" + std::to_string(count);
      }
    };
    add_detail("invalid_gid", invalid_gid_errors);
    add_detail("wrong_rank", wrong_rank_errors);
    add_detail("invalid_species", invalid_species_errors);
    add_detail("nonfinite", nonfinite_errors);
    add_detail("out_of_bounds", bounds_errors);
    add_detail("local_duplicate_tag", local_duplicate_errors);
    add_detail("count_mismatch", count_errors);
    add_detail("global_duplicate_tag", global_duplicate_errors);
    if (!first_local_error.empty()) {
      std::cout << "Particle consistency detail: " << first_local_error << std::endl;
    }
    FatalConsistencyError(label, details);
  }
}

//----------------------------------------------------------------------------------------
// ValidateRelativisticRestartState()
// Enforce the typed-v2 semantic boundary independently of optional debug checks.

void Particles::ValidateRelativisticRestartState(const std::string &label) {
  Mesh *pm = pmy_pack->pmesh;
  auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  auto pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  std::vector<int64_t> identity_keys;
  identity_keys.reserve(std::max(0, nprtcl_thispack));
  RegionSize rs;
  for (int p=0; p<nprtcl_thispack; ++p) {
    int gid = pi_h(PGID,p);
    int tag = pi_h(PTAG,p);
    int spec = pi_h(PSP,p);
    if (gid < 0 || gid >= pm->nmb_total ||
        pm->rank_eachmb[gid] != global_variable::my_rank) {
      FatalConsistencyError(label, "typed-v2 particle has invalid rank ownership");
    }
    if (tag < 0 || spec < 0 || spec >= nspecies) {
      FatalConsistencyError(label, "typed-v2 particle has invalid identity");
    }
    for (int n=0; n<nrdata; ++n) {
      if (!std::isfinite(pr_h(n,p))) {
        FatalConsistencyError(label, "typed-v2 particle has non-finite real state");
      }
    }
    if (pr_h(IPALPHA,p) != alpha_s) {
      FatalConsistencyError(label, "typed-v2 particle signed alpha_s is inconsistent");
    }
    LogicalLocationToRegion(pm, pm->lloc_eachmb[gid], rs);
    if (!ContainsParticle(pm, rs, pm->lloc_eachmb[gid],
                          pr_h(IPX,p), pr_h(IPY,p), pr_h(IPZ,p))) {
      FatalConsistencyError(label, "typed-v2 particle lies outside its saved MeshBlock");
    }
    identity_keys.push_back((static_cast<int64_t>(spec) << 32) ^
                            static_cast<unsigned int>(tag));
  }
  std::sort(identity_keys.begin(), identity_keys.end());
  for (std::size_t n=1; n<identity_keys.size(); ++n) {
    if (identity_keys[n] == identity_keys[n-1]) {
      FatalConsistencyError(label, "typed-v2 particle identity is duplicated");
    }
  }
}

//----------------------------------------------------------------------------------------
// CheckMotionBounds()
// Optional guardrail for the neighbor-table particle exchange path.  The current
// exchange assumes a pushed particle can be assigned through the immediate neighbor
// stencil of its pre-push MeshBlock.

void Particles::CheckMotionBounds(const std::string &label) {
  if (!check_motion_bounds) {return;}

  Mesh *pm = pmy_pack->pmesh;
  auto pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  RegionSize rs;
  for (int p=0; p<nprtcl_thispack; ++p) {
    int gid = pi_h(PGID,p);
    if (gid < 0 || gid >= pm->nmb_total) {
      FatalMotionError(label, "particle " + std::to_string(p) +
                       " has invalid PGID " + std::to_string(gid));
    }
    LogicalLocationToRegion(pm, pm->lloc_eachmb[gid], rs);
    Real lx1 = rs.x1max - rs.x1min;
    Real lx2 = rs.x2max - rs.x2min;
    Real lx3 = rs.x3max - rs.x3min;
    bool bad_x1 = (pr_h(IPX,p) < rs.x1min - lx1) ||
                  (pr_h(IPX,p) > rs.x1max + lx1);
    bool bad_x2 = pm->multi_d && ((pr_h(IPY,p) < rs.x2min - lx2) ||
                                  (pr_h(IPY,p) > rs.x2max + lx2));
    bool bad_x3 = pm->three_d && ((pr_h(IPZ,p) < rs.x3min - lx3) ||
                                  (pr_h(IPZ,p) > rs.x3max + lx3));
    if (bad_x1 || bad_x2 || bad_x3) {
      std::string msg = "species=" + std::to_string(pi_h(PSP,p)) +
                        " tag=" + std::to_string(pi_h(PTAG,p)) +
                        " PGID=" + std::to_string(gid) +
                        " position=(" + std::to_string(pr_h(IPX,p)) + ", " +
                        std::to_string(pr_h(IPY,p)) + ", " +
                        std::to_string(pr_h(IPZ,p)) + ") allowed x1=[" +
                        std::to_string(rs.x1min - lx1) + ", " +
                        std::to_string(rs.x1max + lx1) + "]";
      FatalMotionError(label, msg);
    }
  }
}

//----------------------------------------------------------------------------------------
// LogPerformance()
// Optional local performance counter used for CR tracer tests and diagnostics.

void Particles::LogPerformance(const std::string &label, int64_t npushed,
                               int64_t nsent, int64_t nrecv, int64_t nremapped,
                               int64_t ndiag, int64_t nmessages, int64_t nbytes) {
  if (!log_performance) {return;}

  int64_t local[7] = {npushed, nsent, nrecv, nremapped, ndiag, nmessages, nbytes};
  int64_t global[7] = {0, 0, 0, 0, 0, 0, 0};
#if MPI_PARALLEL_ENABLED
  MPI_Reduce(local, global, 7, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  for (int n=0; n<7; ++n) {global[n] = local[n];}
#endif
  if (global_variable::my_rank == 0) {
    std::cout << "Particle performance " << label
              << ": pushed=" << global[0]
              << " sent=" << global[1]
              << " received=" << global[2]
              << " remapped=" << global[3]
              << " diagnostics=" << global[4]
              << " messages=" << global[5]
              << " bytes=" << global[6] << std::endl;
  }
}

//----------------------------------------------------------------------------------------
// RemapAfterAMR()
// AMR changes the leaf MeshBlock list and may move blocks between ranks.  Particles are
// remapped by position onto the new leaf list, then the existing particle MPI path moves
// off-rank particles.  This is intentionally separate from face-centered AMR repair.

void Particles::RemapAfterAMR() {
  Mesh *pm = pmy_pack->pmesh;
  ParticlesBoundaryValues *pb = pbval_part;
  int send_capacity = std::max(1, nprtcl_thispack);
  Kokkos::realloc(pb->sendlist, send_capacity);
  pb->nprtcl_send = 0;
  pb->nprtcl_recv = 0;
  int nremapped = nprtcl_thispack;

  if (nprtcl_thispack > 0 &&
      amr_remap_mode == ParticlesAMRRemapMode::device_table) {
    DvceArray1D<int> gid_table;
    DvceArray1D<int> rank_table;
    int nloc1, nloc2, nloc3;
    bool built_table = BuildAMRLookupTable(pm, amr_lookup_max_cells, gid_table,
                                           rank_table, nloc1, nloc2, nloc3);
    if (built_table) {
      auto &pr = prtcl_rdata;
      auto &pi = prtcl_idata;
      auto &slist = pb->sendlist;
      DvceArray1D<int> atom_count("amr_particle_remap_count", 2);
      Kokkos::deep_copy(atom_count, 0);
      int myrank = global_variable::my_rank;
      bool multi_d = pm->multi_d;
      bool three_d = pm->three_d;
      Real x1min = pm->mesh_size.x1min;
      Real x1max = pm->mesh_size.x1max;
      Real x2min = pm->mesh_size.x2min;
      Real x2max = pm->mesh_size.x2max;
      Real x3min = pm->mesh_size.x3min;
      Real x3max = pm->mesh_size.x3max;

      par_for("prtcl_amr_remap_table",DevExeSpace(),0,(nprtcl_thispack-1),
      KOKKOS_LAMBDA(const int p) {
        Real x1 = pr(IPX,p);
        Real x2 = pr(IPY,p);
        Real x3 = pr(IPZ,p);
        WrapCoordinateDevice(x1, x1min, x1max);
        if (multi_d) {WrapCoordinateDevice(x2, x2min, x2max);}
        if (three_d) {WrapCoordinateDevice(x3, x3min, x3max);}

        int ix = LocationIndexDevice(x1, x1min, x1max, nloc1);
        int iy = multi_d ? LocationIndexDevice(x2, x2min, x2max, nloc2) : 0;
        int iz = three_d ? LocationIndexDevice(x3, x3min, x3max, nloc3) : 0;
        std::int64_t idx = static_cast<std::int64_t>(ix) +
                           static_cast<std::int64_t>(nloc1)*
                           (static_cast<std::int64_t>(iy) +
                            static_cast<std::int64_t>(nloc2)*
                            static_cast<std::int64_t>(iz));
        int dest_gid = gid_table(idx);
        int dest_rank = rank_table(idx);
        if (dest_gid < 0 || dest_rank < 0) {
          Kokkos::atomic_fetch_add(&atom_count(1), 1);
          return;
        }

        pr(IPX,p) = x1;
        pr(IPY,p) = x2;
        pr(IPZ,p) = x3;
        pi(PGID,p) = dest_gid;
#if MPI_PARALLEL_ENABLED
        if (dest_rank != myrank) {
          int index = Kokkos::atomic_fetch_add(&atom_count(0), 1);
          slist.d_view(index).prtcl_indx = p;
          slist.d_view(index).dest_gid = dest_gid;
          slist.d_view(index).dest_rank = dest_rank;
        }
#endif
      });
      HostArray1D<int> count_h("amr_particle_remap_count_h", 2);
      Kokkos::deep_copy(count_h, atom_count);
      if (count_h(1) != 0) {
        FatalAMRRemapError(std::to_string(count_h(1)) +
                           " particles mapped outside the AMR lookup table");
      }
      pb->nprtcl_send = count_h(0);
      Kokkos::resize(pb->sendlist, pb->nprtcl_send);
      pb->sendlist.template modify<DevExeSpace>();
      pb->sendlist.template sync<HostMemSpace>();
    } else {
      if (global_variable::my_rank == 0) {
        std::cout << "### WARNING: <particles>/amr_remap=device_table requested, "
                  << "but the AMR lookup table would exceed "
                  << amr_lookup_max_cells << " entries; using host_tree"
                  << std::endl;
      }
      amr_remap_mode = ParticlesAMRRemapMode::host_tree;
    }
  }

  if (nprtcl_thispack > 0 &&
      amr_remap_mode == ParticlesAMRRemapMode::host_tree) {
    auto x1_d = Kokkos::subview(prtcl_rdata, static_cast<int>(IPX), Kokkos::ALL());
    auto x2_d = Kokkos::subview(prtcl_rdata, static_cast<int>(IPY), Kokkos::ALL());
    auto x3_d = Kokkos::subview(prtcl_rdata, static_cast<int>(IPZ), Kokkos::ALL());
    auto gid_d = Kokkos::subview(prtcl_idata, static_cast<int>(PGID), Kokkos::ALL());
    HostArray1D<Real> x1_h("prtcl_remap_x1", nprtcl_thispack);
    HostArray1D<Real> x2_h("prtcl_remap_x2", nprtcl_thispack);
    HostArray1D<Real> x3_h("prtcl_remap_x3", nprtcl_thispack);
    HostArray1D<int> gid_h("prtcl_remap_gid", nprtcl_thispack);
    Kokkos::deep_copy(x1_h, x1_d);
    Kokkos::deep_copy(x2_h, x2_d);
    Kokkos::deep_copy(x3_h, x3_d);
    Kokkos::deep_copy(gid_h, gid_d);

    int myrank = global_variable::my_rank;
    for (int p=0; p<nprtcl_thispack; ++p) {
      WrapCoordinate(x1_h(p), pm->mesh_size.x1min, pm->mesh_size.x1max);
      if (pm->multi_d) {
        WrapCoordinate(x2_h(p), pm->mesh_size.x2min, pm->mesh_size.x2max);
      }
      if (pm->three_d) {
        WrapCoordinate(x3_h(p), pm->mesh_size.x3min, pm->mesh_size.x3max);
      }

      int dest_gid = pm->FindMeshBlockContainingPosition(x1_h(p), x2_h(p), x3_h(p));
      if (dest_gid < 0) {
        FatalAMRRemapError("could not map particle " + std::to_string(p) +
                           " after AMR");
      }
      gid_h(p) = dest_gid;
#if MPI_PARALLEL_ENABLED
      int dest_rank = pm->rank_eachmb[dest_gid];
      if (dest_rank != myrank) {
        int index = pb->nprtcl_send++;
        pb->sendlist.h_view(index).prtcl_indx = p;
        pb->sendlist.h_view(index).dest_gid = dest_gid;
        pb->sendlist.h_view(index).dest_rank = dest_rank;
      }
#endif
    }

    Kokkos::deep_copy(x1_d, x1_h);
    Kokkos::deep_copy(x2_d, x2_h);
    Kokkos::deep_copy(x3_d, x3_h);
    Kokkos::deep_copy(gid_d, gid_h);
    Kokkos::resize(pb->sendlist, pb->nprtcl_send);
    pb->sendlist.template modify<HostMemSpace>();
    pb->sendlist.template sync<DevExeSpace>();
  }

  if (validate_amr_lookup ||
      consistency_mode == ParticlesConsistencyMode::local ||
      consistency_mode == ParticlesConsistencyMode::full) {
    auto pr_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
    auto pi_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
    for (int p=0; p<nprtcl_thispack; ++p) {
      int tree_gid = pm->FindMeshBlockContainingPosition(pr_h(IPX,p),
                                                         pr_h(IPY,p),
                                                         pr_h(IPZ,p));
      int gid = pi_h(PGID,p);
      if (tree_gid != gid) {
        FatalAMRRemapError("lookup validation mismatch for particle " +
                           std::to_string(p) + ": table gid=" +
                           std::to_string(gid) + " tree gid=" +
                           std::to_string(tree_gid));
      }
      RegionSize rs;
      LogicalLocationToRegion(pm, pm->lloc_eachmb[gid], rs);
      if (!ContainsParticle(pm, rs, pm->lloc_eachmb[gid],
                            pr_h(IPX,p), pr_h(IPY,p), pr_h(IPZ,p))) {
        FatalAMRRemapError("lookup validation selected MeshBlock " +
                           std::to_string(gid) + " that does not contain " +
                           "particle " + std::to_string(p));
      }
    }
  }
  if (nprtcl_thispack == 0) {
    Kokkos::resize(pb->sendlist, 0);
    pb->sendlist.template modify<HostMemSpace>();
    pb->sendlist.template sync<DevExeSpace>();
  }

#if MPI_PARALLEL_ENABLED
  pb->CountSendsAndRecvs();
  pb->InitPrtclRecv();
  pb->PackAndSendPrtcls();
  while (pb->RecvAndUnpackPrtcls() == TaskStatus::incomplete) {}
  pb->ClearPrtclRecv();
  pb->ClearPrtclSend();
#else
  pm->UpdateParticleCounts();
#endif
  LogPerformance("amr_remap", 0, pb->nprtcl_send, pb->nprtcl_recv, nremapped, 0);
  CheckConsistency("AMR remap");
  MarkRelativisticTimestepBoundDirty();
}

} // namespace particles
