//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.cpp
//! \brief implementation of Particles class constructor and assorted other functions

#include <iostream>
#include <string>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "bvals/bvals.hpp"
#include "particles.hpp"
#include "srcterms/srcterms.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace particles {
namespace {

std::uint64_t GetOrAddUInt64(ParameterInput *pin, const std::string &block,
                             const std::string &name, std::uint64_t value) {
  std::string text = pin->GetOrAddString(block, name, std::to_string(value));
  if (text.empty() || text.front() == '-') {
    std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
              << block << "/" << name << " must be an unsigned 64-bit integer"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::size_t consumed = 0;
  try {
    std::uint64_t parsed = std::stoull(text, &consumed, 10);
    if (consumed != text.size()) throw std::invalid_argument("trailing characters");
    return parsed;
  } catch (const std::exception &) {
    std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
              << block << "/" << name << "='" << text
              << "' is not an unsigned 64-bit integer" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

} // namespace

//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Particles::Particles(MeshBlockPack *ppack, ParameterInput *pin,
                     bool infer_covariance_model_from_restart) :
    pmy_pack(ppack) {
  // check this is at least a 2D problem
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle module only works in 2D/3D" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // select particle type
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("cosmic_ray") == 0) {
      particle_type = ParticleType::cosmic_ray;
    } else if (ptype.compare("lagrangian_mc") == 0) {
      particle_type = ParticleType::lagrangian_mc;
    } else if (ptype.compare("lagrangian_ito") == 0) {
      particle_type = ParticleType::lagrangian_ito;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle type = '" << ptype << "' not recognized"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // select pusher algorithm
  {
    std::string ppush = pin->GetString("particles","pusher");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
    } else if (ppush.compare("lagrangian_mc") == 0) {
      pusher = ParticlesPusher::lagrangian_mc;
    } else if (ppush.compare("ito2") == 0) {
      pusher = ParticlesPusher::ito2;
    } else if (ppush.compare("ito3") == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "pusher=ito3 is intentionally not implemented"
                << std::endl;
      std::exit(EXIT_FAILURE);
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher must be specified in <particles> block"
                <<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (particle_type == ParticleType::lagrangian_mc &&
      pusher != ParticlesPusher::lagrangian_mc) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "particle_type=lagrangian_mc requires "
              << "pusher=lagrangian_mc" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (particle_type == ParticleType::lagrangian_ito &&
      pusher != ParticlesPusher::ito2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "particle_type=lagrangian_ito requires "
              << "pusher=ito2; Ito-3 is intentionally not implemented" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (particle_type != ParticleType::lagrangian_ito &&
      pusher == ParticlesPusher::ito2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "pusher=ito2 requires particle_type=lagrangian_ito"
              << std::endl;
    std::exit(EXIT_FAILURE);
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
        Real ppc = pin->GetOrAddReal("particles","ppc",1.0);
        auto &indcs = pmy_pack->pmesh->mb_indcs;
        int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
        Real r_npart = ppc*static_cast<Real>((pmy_pack->nmb_thispack)*ncells);
        nprtcl_thispack = static_cast<int>(r_npart);
        int ndim=4;
        if (pmy_pack->pmesh->three_d) {ndim+=2;}
        nrdata = ndim;
        nidata = 2;
        break;
      }
    case ParticleType::lagrangian_mc:
    case ParticleType::lagrangian_ito:
      {
        if (pmy_pack->phydro == nullptr && pmy_pack->pmhd == nullptr) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "flux tracer particles require hydro or mhd"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        std::string fluid_block = (pmy_pack->phydro != nullptr) ? "hydro" : "mhd";
        std::string source_block = fluid_block + "_srcterms";
        if (SourceTerms::MassChangeDeclared(source_block, pin)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "flux tracer particles do not support source terms "
                    << "that change mass density" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        EquationOfState *peos = (pmy_pack->phydro != nullptr) ?
                                pmy_pack->phydro->peos : pmy_pack->pmhd->peos;
        if (particle_type == ParticleType::lagrangian_mc && !peos->eos_data.is_ideal) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "lagrangian_mc thermo tracers require an ideal-gas "
                    << "EOS in this implementation" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (particle_type == ParticleType::lagrangian_ito) {
          if (pin->DoesParameterExist("particles", "ito_order") &&
              pin->GetInteger("particles", "ito_order") != 2) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "lagrangian_ito supports only ito_order=2"
                      << std::endl;
            std::exit(EXIT_FAILURE);
          }
          if (pin->DoesParameterExist("particles", "tracer_kick_pdf") &&
              pin->GetString("particles", "tracer_kick_pdf") != "uniform") {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "lagrangian_ito uses only the bounded uniform "
                      << "Ito-2 kick distribution" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          std::string covariance_model = pin->GetOrAddString(
              "particles", "ito_covariance_model", "published_diagonal");
          if (infer_covariance_model_from_restart) {
            infer_ito_covariance_model_from_restart = true;
            ito_ncoeff = ITO_NFULL_COEFF;
          } else if (covariance_model == "published_diagonal") {
            ito_covariance_model = ItoCovarianceModel::published_diagonal;
            ito_ncoeff = ITO_NDIAG_COEFF;
          } else if (covariance_model == "full_finite_step") {
            ito_covariance_model = ItoCovarianceModel::full_finite_step;
            ito_ncoeff = ITO_NFULL_COEFF;
          } else {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line "
                      << __LINE__ << std::endl
                      << "particles/ito_covariance_model must be "
                      << "published_diagonal or full_finite_step" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          std::string evolution = pin->GetString("time", "evolution");
          std::string integrator = pin->GetOrAddString("time", "integrator", "rk2");
          if (evolution != "dynamic" ||
              (integrator != "rk1" && integrator != "rk2" && integrator != "rk3")) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "lagrangian_ito requires dynamic evolution with "
                      << "integrator=rk1, rk2, or rk3" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          if (!pmy_pack->pmesh->strictly_periodic) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "lagrangian_ito currently requires periodic "
                      << "boundaries in every active dimension" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          if (pmy_pack->pcoord->is_special_relativistic ||
              pmy_pack->pcoord->is_general_relativistic ||
              pmy_pack->pcoord->is_dynamical_relativistic) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "lagrangian_ito currently supports only "
                      << "non-relativistic Cartesian fluid evolution" << std::endl;
            std::exit(EXIT_FAILURE);
          }
        }
        nprtcl_thispack = 0;
        nrdata = LMC_NREAL;
        nidata = PSEEDID + 1;
        random_seed = pin->GetOrAddInteger("particles","random_seed",12345);
        next_tracer_tag = GetOrAddUInt64(pin, "particles", "next_tracer_tag", 0);
        ito_probability_target =
            pin->GetOrAddReal("particles", "ito_probability_target", 0.99);
        if (!(ito_probability_target > 0.0 && ito_probability_target <= 1.0)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "particles/ito_probability_target must be in (0,1]"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        ParseTracerSeedSchedules(pin);
        if (pmy_pack->phydro != nullptr) {
          pmy_pack->phydro->SetSaveUFlxIdn();
        } else {
          pmy_pack->pmhd->SetSaveUFlxIdn();
        }
        dtnew = std::numeric_limits<Real>::max();
        if (particle_type == ParticleType::lagrangian_ito) {
          int nmb = std::max((pmy_pack->nmb_thispack),
                             (pmy_pack->pmesh->nmb_maxperrank));
          auto &indcs = pmy_pack->pmesh->mb_indcs;
          int ncells1 = indcs.nx1 + 2*indcs.ng;
          int ncells2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*indcs.ng : 1;
          int ncells3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*indcs.ng : 1;
          Kokkos::realloc(ito_coeff, nmb, ito_ncoeff, ncells3, ncells2, ncells1);
          if (pmy_pack->pmesh->multilevel) {
            int nccells1 = indcs.cnx1 + 2*indcs.ng;
            int nccells2 = (indcs.cnx2 > 1) ? indcs.cnx2 + 2*indcs.ng : 1;
            int nccells3 = (indcs.cnx3 > 1) ? indcs.cnx3 + 2*indcs.ng : 1;
            Kokkos::realloc(coarse_ito_coeff, nmb, ito_ncoeff,
                            nccells3, nccells2, nccells1);
          }
          Kokkos::realloc(ito_invalid, 1);
        }
        break;
      }
    default:
      break;
  }
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);
  Kokkos::realloc(prtcl_tag, nprtcl_thispack);

  // allocate boundary object
  pbval_part = new ParticlesBoundaryValues(this, pin);
  if (IsIto2()) {
    pbval_ito = new MeshBoundaryValuesCC(ppack, pin, false);
    pbval_ito->InitializeBuffers(ito_ncoeff);
  }
}

//----------------------------------------------------------------------------------------
// destructor

Particles::~Particles() {
  if (pbval_ito != nullptr) delete pbval_ito;
}

//----------------------------------------------------------------------------------------
// GetLagrangianMCScalarCount()

int Particles::GetLagrangianMCScalarCount() const {
  if (pmy_pack->phydro != nullptr) return pmy_pack->phydro->nscalars;
  if (pmy_pack->pmhd != nullptr) return pmy_pack->pmhd->nscalars;
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn void Particles::CheckMassFloorCompatibility
//! \brief Fail before a flux-tracer update if the EOS injected gas mass through a floor.

void Particles::CheckMassFloorCompatibility() {
  int floor_count = pmy_pack->pmesh->ecounter.neos_dfloor;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &floor_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (floor_count > 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "an EOS density floor or MHD magnetization ceiling "
              << "injected gas mass before a flux-tracer update; tracer creation or "
              << "weight adjustment for that mass is not implemented" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
// CreateParticleTags()
// Assigns tags to particles (unique integer).  Note that tracked particles are always
// those with tag numbers less than ntrack.

void Particles::CreateParticleTags(ParameterInput *pin) {
  if (IsFluxTracer()) return;

  std::string assign = pin->GetOrAddString("particles","assign_tag","index_order");

  // tags are assigned sequentially within this rank, starting at 0 with rank=0
  if (assign.compare("index_order") == 0) {
    std::uint64_t tagstart = 0;
    for (int n=1; n<=global_variable::my_rank; ++n) {
      tagstart += static_cast<std::uint64_t>(
          pmy_pack->pmesh->nprtcl_eachrank[n-1]);
    }

    auto &ptag = prtcl_tag;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      ptag(p) = tagstart + static_cast<std::uint64_t>(p);
    });

  // tags are assigned sequentially across ranks
  } else if (assign.compare("rank_order") == 0) {
    int myrank = global_variable::my_rank;
    int nranks = global_variable::nranks;
    auto &ptag = prtcl_tag;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      ptag(p) = static_cast<std::uint64_t>(myrank) +
                 static_cast<std::uint64_t>(nranks)*static_cast<std::uint64_t>(p);
    });

  // tag algorithm not recognized, so quit with error
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle tag assignment type = '" << assign << "' not recognized"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

} // namespace particles
