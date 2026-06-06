//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pgen.cpp
//! \brief Implementation of constructors and functions in class ProblemGenerator.
//! Default constructor calls problem generator function, while  constructor for restarts
//! reads data from restart file, as well as re-initializing problem-specific data.

#include <array>
#include <cstdio>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <algorithm>
#include <vector>

#include "athena.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "srcterms/turb_driver.hpp"
#include "particles/particles.hpp"
#include "outputs/restart_utils.hpp"
#include "pgen.hpp"

namespace {

constexpr IOWrapperSizeT kMaxSupportedRestartParticlePayloadBytesPerRank =
    static_cast<IOWrapperSizeT>(8) << 30;

struct RestartBlockRequest {
  int local_index;
  int global_id;
};

struct ParticleRestartSectionMeta {
  int version = -1;
  int nmb_section = 0;
  int nrdata = 0;
  int nidata = 0;
  int rst_nout1 = 0;
  int rst_nout2 = 0;
  int rst_nout3 = 0;
  int has_moments = 0;
  int has_edge = 0;
  int moment_cnt = 0;
  int edge1_cnt = 0;
  int edge2_cnt = 0;
  int edge3_cnt = 0;
  int state_kind = -1;
  int physical_mode = -1;
  Real cr_light_speed = 0.0;
  std::array<int, particles::Particles::NPIC_RESTART_MODEL_INTS> model_ints;
  std::array<Real, particles::Particles::NPIC_RESTART_MODEL_REALS> model_reals;
  IOWrapperSizeT npart_section = 0;
  IOWrapperSizeT pr_real_offset = 0;
  IOWrapperSizeT pr_int_offset = 0;
  IOWrapperSizeT moments_offset = 0;
  IOWrapperSizeT edge1_offset = 0;
  IOWrapperSizeT edge2_offset = 0;
  IOWrapperSizeT edge3_offset = 0;
  std::vector<int> mb_counts;
  std::vector<IOWrapperSizeT> mb_offsets;
};

void AbortParticleRestartLayoutOverflow() {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl
            << "Particle restart layout exceeds supported offset or allocation limits."
            << std::endl;
  restart_utils::AbortOnFatalError();
}

IOWrapperSizeT CheckedParticleRestartProduct(const IOWrapperSizeT count,
                                             const IOWrapperSizeT width) {
  constexpr IOWrapperSizeT max = std::numeric_limits<IOWrapperSizeT>::max();
  if (count != 0 && width > max/count) AbortParticleRestartLayoutOverflow();
  return count*width;
}

IOWrapperSizeT CheckedParticleRestartProduct(
    std::initializer_list<IOWrapperSizeT> factors) {
  IOWrapperSizeT product = 1;
  for (const IOWrapperSizeT factor : factors) {
    product = CheckedParticleRestartProduct(product, factor);
  }
  return product;
}

void ValidateParticleRestartAllocation(std::initializer_list<int> extents) {
  std::size_t count = 1;
  constexpr std::size_t max = std::numeric_limits<std::size_t>::max();
  constexpr std::size_t element_bytes =
      sizeof(Real) > sizeof(int) ? sizeof(Real) : sizeof(int);
  for (const int extent : extents) {
    if (extent < 0 || (extent != 0 &&
                       count > max/static_cast<std::size_t>(extent))) {
      AbortParticleRestartLayoutOverflow();
    }
    count *= static_cast<std::size_t>(extent);
  }
  if (count > max/element_bytes) AbortParticleRestartLayoutOverflow();
}

void AdvanceParticleRestartOffset(IOWrapperSizeT &offset, const IOWrapperSizeT count,
                                  const IOWrapperSizeT width) {
  constexpr IOWrapperSizeT max = std::numeric_limits<IOWrapperSizeT>::max();
  const IOWrapperSizeT increment = CheckedParticleRestartProduct(count, width);
  if (offset > max - increment) AbortParticleRestartLayoutOverflow();
  offset += increment;
}

void AdvanceParticleRestartRealArrayOffset(
    IOWrapperSizeT &offset, std::initializer_list<int> extents) {
  IOWrapperSizeT count = 1;
  for (const int extent : extents) {
    if (extent < 0) AbortParticleRestartLayoutOverflow();
    count = CheckedParticleRestartProduct(count, static_cast<IOWrapperSizeT>(extent));
  }
  AdvanceParticleRestartOffset(offset, count, sizeof(Real));
}

void AdvanceParticleRestartPastFaceArrayRemainder(
    IOWrapperSizeT &offset, const IOWrapperSizeT data_size,
    std::initializer_list<IOWrapperSizeT> face_counts) {
  IOWrapperSizeT face_bytes = 0;
  for (const IOWrapperSizeT count : face_counts) {
    AdvanceParticleRestartOffset(face_bytes, count, sizeof(Real));
  }
  if (face_bytes > data_size) AbortParticleRestartLayoutOverflow();
  AdvanceParticleRestartOffset(offset, data_size - face_bytes, 1);
}

IOWrapperSizeT CheckedParticleRestartDataOffset(const IOWrapperSizeT base,
                                                const IOWrapperSizeT relative) {
  IOWrapperSizeT offset = base;
  AdvanceParticleRestartOffset(offset, relative, 1);
  return offset;
}

int CheckedParticleRestartExtentPlusOne(const int extent) {
  if (extent < 0 || extent == std::numeric_limits<int>::max()) {
    AbortParticleRestartLayoutOverflow();
  }
  return extent + 1;
}

int CheckedParticleRestartOutputExtent(const int active_cells, const int ghost_cells) {
  if (active_cells <= 0 || ghost_cells < 0) AbortParticleRestartLayoutOverflow();
  IOWrapperSizeT extent = static_cast<IOWrapperSizeT>(active_cells);
  if (active_cells > 1) {
    AdvanceParticleRestartOffset(extent, static_cast<IOWrapperSizeT>(ghost_cells), 2);
  }
  if (extent > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max())) {
    AbortParticleRestartLayoutOverflow();
  }
  return static_cast<int>(extent);
}

int CheckedParticleRestartIntProduct(std::initializer_list<int> factors) {
  IOWrapperSizeT product = 1;
  for (const int factor : factors) {
    if (factor < 0) AbortParticleRestartLayoutOverflow();
    product = CheckedParticleRestartProduct(product, static_cast<IOWrapperSizeT>(factor));
  }
  if (product > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max())) {
    AbortParticleRestartLayoutOverflow();
  }
  return static_cast<int>(product);
}

int CheckedLocalParticleCount(const IOWrapperSizeT count, const int nrdata,
                              const int nidata) {
  if (nrdata <= 0 || nidata <= 0) AbortParticleRestartLayoutOverflow();
  const IOWrapperSizeT real_values =
      CheckedParticleRestartProduct(count, static_cast<IOWrapperSizeT>(nrdata));
  const IOWrapperSizeT int_values =
      CheckedParticleRestartProduct(count, static_cast<IOWrapperSizeT>(nidata));
  const IOWrapperSizeT int_bytes = CheckedParticleRestartProduct(int_values, sizeof(int));
  const IOWrapperSizeT real_bytes =
      CheckedParticleRestartProduct(real_values, sizeof(Real));
  if (count > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max()) ||
      real_values > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max()) ||
      int_bytes > static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max()) ||
      real_bytes > kMaxSupportedRestartParticlePayloadBytesPerRank - int_bytes) {
    AbortParticleRestartLayoutOverflow();
  }
  return static_cast<int>(count);
}

template <typename IntArray>
void ValidateRestoredParticleIDData(const particles::Particles *ppart,
                                    const IntArray &h_pi, const int p,
                                    const int gids_local, const int nmb_local) {
  const int gid = h_pi(PGID, p);
  if (gid < gids_local || gid >= (gids_local + nmb_local)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restarted particle gid is not local after restore (gid="
              << gid << ")." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (ppart->particle_type == ParticleType::cosmic_ray) {
    const int sp = h_pi(PSP, p);
    if (sp < 0 || sp >= ppart->nspecies) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restarted cosmic-ray particle species is out of range "
                << "(species=" << sp << ", nspecies=" << ppart->nspecies << ")."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
}

void LoadParticleRestartDataSingleFile(Mesh *pm,
                                       IOWrapperSizeT headeroffset,
                                       IOWrapperSizeT data_stride,
                                       int nout1, int nout2, int nout3) {
  auto *ppart = pm->pmb_pack->ppart;
  if (ppart == nullptr) return;

  constexpr std::uint64_t kPicMagic = 0x5049435253543031ULL;
  constexpr int kPicVersion = particles::Particles::PIC_RESTART_SCHEMA_VERSION;

  const RestartMetaData &meta = pm->restart_meta;
  if (meta.file_name.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart metadata missing file name for single-file restart."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (meta.rank_eachmb.size() != static_cast<std::size_t>(pm->nmb_total)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart metadata inconsistent with MeshBlock count."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (meta.original_nranks <= 0 ||
      meta.gids_eachrank.size() != static_cast<std::size_t>(meta.original_nranks) ||
      meta.nmb_eachrank.size() != static_cast<std::size_t>(meta.original_nranks)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart metadata missing original rank layout."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  MeshBlockPack *pack = pm->pmb_pack;
  const int nmb_local = pack->nmb_thispack;
  const int nrdata = ppart->nrdata;
  const int nidata = ppart->nidata;
  const int expected_moment_cnt = CheckedParticleRestartIntProduct(
      {particles::Particles::NMOM, nout3, nout2, nout1});
  const int nout1p1 = CheckedParticleRestartExtentPlusOne(nout1);
  const int nout2p1 = CheckedParticleRestartExtentPlusOne(nout2);
  const int nout3p1 = CheckedParticleRestartExtentPlusOne(nout3);
  const int expected_edge1_cnt =
      CheckedParticleRestartIntProduct({nout3p1, nout2p1, nout1});
  const int expected_edge2_cnt =
      CheckedParticleRestartIntProduct({nout3p1, nout2, nout1p1});
  const int expected_edge3_cnt =
      CheckedParticleRestartIntProduct({nout3, nout2p1, nout1p1});

  std::vector<std::vector<RestartBlockRequest>> requests(meta.original_nranks);
  for (int m=0; m<nmb_local; ++m) {
    int gid = pack->pmb->mb_gid.h_view(m);
    if (gid < 0 || gid >= pm->nmb_total) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Invalid MeshBlock gid encountered during restart."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    int src_rank = meta.rank_eachmb[gid];
    if (src_rank < 0 || src_rank >= meta.original_nranks) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart metadata contains invalid rank assignments."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    requests[src_rank].push_back({m, gid});
  }

  std::vector<std::string> rank_paths(meta.original_nranks);
  const std::string base_prefix = meta.base_dir.empty() ? "" :
      (meta.base_dir == "/" ? "/" : meta.base_dir + "/");
  for (int r=0; r<meta.original_nranks; ++r) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", r);
    rank_paths[r] = base_prefix + rank_dir + "/" + meta.file_name;
  }

  std::vector<int> local_mb_counts(nmb_local, 0);
  std::vector<ParticleRestartSectionMeta> src_meta(meta.original_nranks);
  std::vector<bool> src_used(meta.original_nranks, false);
  int ref_has_moments = -1;
  int ref_has_edge = -1;
  int ref_moment_cnt = 0;
  int ref_edge1_cnt = 0;
  int ref_edge2_cnt = 0;
  int ref_edge3_cnt = 0;
  int ref_state_kind = -1;
  int ref_physical_mode = -1;
  Real ref_cr_light_speed = 0.0;
  std::array<int, particles::Particles::NPIC_RESTART_MODEL_INTS> ref_model_ints;
  std::array<Real, particles::Particles::NPIC_RESTART_MODEL_REALS> ref_model_reals;

  for (int r=0; r<meta.original_nranks; ++r) {
    auto &reqs = requests[r];

    IOWrapper srcfile;
    srcfile.Open(rank_paths[r].c_str(), IOWrapper::FileMode::read, true);

    IOWrapperSizeT section_offset = headeroffset;
    AdvanceParticleRestartOffset(section_offset,
                                 static_cast<IOWrapperSizeT>(meta.nmb_eachrank[r]),
                                 data_stride);

    std::uint64_t pic_magic = 0;
    if (srcfile.Read_bytes_at(&pic_magic, 1, sizeof(std::uint64_t), section_offset,
                              true) != sizeof(std::uint64_t)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
              << "Particle restart state is missing from source restart "
                << "file '" << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (pic_magic != kPicMagic) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Source restart file '" << rank_paths[r]
                << "' has an unknown particle restart marker."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }

    ParticleRestartSectionMeta sm;
    IOWrapperSizeT rd_offset = section_offset;
    AdvanceParticleRestartOffset(rd_offset, 1, sizeof(std::uint64_t));
    auto read_int_meta = [&](int &val, const char *name) {
      if (srcfile.Read_bytes_at(&val, sizeof(int), 1, rd_offset, true) != 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to read particle restart metadata field '" << name
                  << "' from source restart file '" << rank_paths[r] << "'."
                  << std::endl;
        restart_utils::AbortOnFatalError();
      }
      AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
    };

    read_int_meta(sm.version, "version");
    if (sm.version != kPicVersion) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Unsupported particle restart version " << sm.version
                << " in source restart file '" << rank_paths[r] << "'."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    read_int_meta(sm.nmb_section, "nmb_section");
    read_int_meta(sm.nrdata, "nrdata");
    read_int_meta(sm.nidata, "nidata");
    if (sm.nrdata <= 0 || sm.nidata <= 0) AbortParticleRestartLayoutOverflow();
    read_int_meta(sm.rst_nout1, "nout1");
    read_int_meta(sm.rst_nout2, "nout2");
    read_int_meta(sm.rst_nout3, "nout3");
    read_int_meta(sm.has_moments, "has_moments");
    read_int_meta(sm.has_edge, "has_edge");
    read_int_meta(sm.moment_cnt, "moment_cnt");
    read_int_meta(sm.edge1_cnt, "edge1_cnt");
    read_int_meta(sm.edge2_cnt, "edge2_cnt");
    read_int_meta(sm.edge3_cnt, "edge3_cnt");
    read_int_meta(sm.state_kind, "state_kind");
    read_int_meta(sm.physical_mode, "physical_mode");
    if (srcfile.Read_Reals_at(&sm.cr_light_speed, 1, rd_offset, true) != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart metadata field 'cr_light_speed' "
                << "from source restart file '" << rank_paths[r] << "'."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    AdvanceParticleRestartOffset(rd_offset, 1, sizeof(Real));
    if (srcfile.Read_bytes_at(sm.model_ints.data(), sizeof(int), sm.model_ints.size(),
                              rd_offset, true) != sm.model_ints.size()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart extension-model integer metadata "
                << "from source restart file '" << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    AdvanceParticleRestartOffset(rd_offset, sm.model_ints.size(), sizeof(int));
    if (srcfile.Read_Reals_at(sm.model_reals.data(), sm.model_reals.size(), rd_offset,
                              true) != sm.model_reals.size()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart extension-model numerical metadata "
                << "from source restart file '" << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    AdvanceParticleRestartOffset(rd_offset, sm.model_reals.size(), sizeof(Real));

    if (srcfile.Read_bytes_at(&sm.npart_section, sizeof(IOWrapperSizeT), 1, rd_offset,
                              true) != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle-count metadata from source "
                << "restart file '" << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    AdvanceParticleRestartOffset(rd_offset, 1, sizeof(IOWrapperSizeT));

    if (sm.nmb_section != meta.nmb_eachrank[r]) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart MeshBlock count mismatch in source restart file '"
                << rank_paths[r] << "' (file=" << sm.nmb_section << ", metadata="
                << meta.nmb_eachrank[r] << ")." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (sm.nrdata != nrdata || sm.nidata != nidata) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart data layout mismatch in source restart "
                << "file '" << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    const int expected_state_kind = ppart->UsesRelativisticCRState() ? 1 : 0;
    const int expected_physical_mode = static_cast<int>(ppart->pic_physical_mode);
    if (sm.state_kind != expected_state_kind ||
        sm.physical_mode != expected_physical_mode ||
        sm.cr_light_speed != ppart->pic_cr_light_speed ||
        !ppart->MatchesRestartModelMetadata(sm.model_ints, sm.model_reals) ||
        !ppart->RestoreRestartModelState(sm.model_reals)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart physical-model metadata mismatch in source "
                << "restart file '" << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (sm.rst_nout1 != nout1 || sm.rst_nout2 != nout2 || sm.rst_nout3 != nout3) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart mesh extents mismatch in source restart file '"
                << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if ((sm.has_moments != 0 && sm.has_moments != 1) ||
        (sm.has_edge != 0 && sm.has_edge != 1)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart optional-array flags are invalid in source "
                << "restart file '" << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (sm.moment_cnt != expected_moment_cnt) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart moment size mismatch in source restart file '"
                << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (sm.has_edge != 0 &&
        (sm.edge1_cnt != expected_edge1_cnt || sm.edge2_cnt != expected_edge2_cnt ||
         sm.edge3_cnt != expected_edge3_cnt)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart edge-current size mismatch in source restart file '"
                << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }

    if (ref_has_moments < 0) {
      ref_has_moments = sm.has_moments;
      ref_has_edge = sm.has_edge;
      ref_moment_cnt = sm.moment_cnt;
      ref_edge1_cnt = sm.edge1_cnt;
      ref_edge2_cnt = sm.edge2_cnt;
      ref_edge3_cnt = sm.edge3_cnt;
      ref_state_kind = sm.state_kind;
      ref_physical_mode = sm.physical_mode;
      ref_cr_light_speed = sm.cr_light_speed;
      ref_model_ints = sm.model_ints;
      ref_model_reals = sm.model_reals;
    } else if (sm.has_moments != ref_has_moments ||
               sm.has_edge != ref_has_edge ||
               sm.moment_cnt != ref_moment_cnt ||
               sm.edge1_cnt != ref_edge1_cnt ||
               sm.edge2_cnt != ref_edge2_cnt ||
               sm.edge3_cnt != ref_edge3_cnt ||
               sm.state_kind != ref_state_kind ||
               sm.physical_mode != ref_physical_mode ||
               sm.cr_light_speed != ref_cr_light_speed ||
               sm.model_ints != ref_model_ints ||
               sm.model_reals != ref_model_reals) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart metadata is inconsistent across source restart "
                << "files." << std::endl;
      restart_utils::AbortOnFatalError();
    }

    sm.mb_counts.assign(sm.nmb_section, 0);
    if (sm.nmb_section > 0) {
      if (srcfile.Read_bytes_at(sm.mb_counts.data(), sizeof(int), sm.nmb_section,
                                rd_offset, true) !=
          static_cast<std::size_t>(sm.nmb_section)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to read particle restart MeshBlock counts from source "
                  << "restart file '" << rank_paths[r] << "'." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
    AdvanceParticleRestartOffset(rd_offset, static_cast<IOWrapperSizeT>(sm.nmb_section),
                                 sizeof(int));

    sm.mb_offsets.assign(sm.nmb_section + 1, 0);
    for (int m=0; m<sm.nmb_section; ++m) {
      if (sm.mb_counts[m] < 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Particle restart MeshBlock count table contains a negative "
                  << "count in source restart file '" << rank_paths[r] << "'."
                  << std::endl;
        restart_utils::AbortOnFatalError();
      }
      sm.mb_offsets[m + 1] = sm.mb_offsets[m];
      AdvanceParticleRestartOffset(sm.mb_offsets[m + 1],
                                   static_cast<IOWrapperSizeT>(sm.mb_counts[m]), 1);
    }
    if (sm.mb_offsets[sm.nmb_section] != sm.npart_section) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart count table is inconsistent in source "
                << "restart file '" << rank_paths[r] << "'." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    (void)CheckedLocalParticleCount(sm.npart_section, sm.nrdata, sm.nidata);

    sm.pr_real_offset = rd_offset;
    AdvanceParticleRestartOffset(rd_offset, sm.npart_section,
                                 CheckedParticleRestartProduct(
                                     static_cast<IOWrapperSizeT>(sm.nrdata),
                                     sizeof(Real)));
    sm.pr_int_offset = rd_offset;
    AdvanceParticleRestartOffset(rd_offset, sm.npart_section,
                                 CheckedParticleRestartProduct(
                                     static_cast<IOWrapperSizeT>(sm.nidata),
                                     sizeof(int)));
    sm.moments_offset = rd_offset;
    if (sm.has_moments != 0) {
      AdvanceParticleRestartOffset(
          rd_offset, static_cast<IOWrapperSizeT>(sm.nmb_section),
          CheckedParticleRestartProduct(
              static_cast<IOWrapperSizeT>(sm.moment_cnt), sizeof(Real)));
    }
    sm.edge1_offset = rd_offset;
    if (sm.has_edge != 0) {
      AdvanceParticleRestartOffset(
          rd_offset, static_cast<IOWrapperSizeT>(sm.nmb_section),
          CheckedParticleRestartProduct(
              static_cast<IOWrapperSizeT>(sm.edge1_cnt), sizeof(Real)));
    }
    sm.edge2_offset = rd_offset;
    if (sm.has_edge != 0) {
      AdvanceParticleRestartOffset(
          rd_offset, static_cast<IOWrapperSizeT>(sm.nmb_section),
          CheckedParticleRestartProduct(
              static_cast<IOWrapperSizeT>(sm.edge2_cnt), sizeof(Real)));
    }
    sm.edge3_offset = rd_offset;
    if (sm.has_edge != 0) {
      AdvanceParticleRestartOffset(
          rd_offset, static_cast<IOWrapperSizeT>(sm.nmb_section),
          CheckedParticleRestartProduct(
              static_cast<IOWrapperSizeT>(sm.edge3_cnt), sizeof(Real)));
    }
    if (rd_offset > srcfile.GetSize(true)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart section exceeds source artifact bounds."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }

    for (const auto &req : reqs) {
      const int src_local = req.global_id - meta.gids_eachrank[r];
      if (src_local < 0 || src_local >= sm.nmb_section) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Restart metadata is inconsistent with particle section "
                  << "layout in source restart file '" << rank_paths[r] << "'."
                  << std::endl;
        restart_utils::AbortOnFatalError();
      }
      local_mb_counts[req.local_index] = sm.mb_counts[src_local];
    }

    src_meta[r] = std::move(sm);
    src_used[r] = true;
    srcfile.Close(true);
  }

  std::vector<IOWrapperSizeT> local_mb_offsets(nmb_local + 1, 0);
  for (int m=0; m<nmb_local; ++m) {
    local_mb_offsets[m + 1] = local_mb_offsets[m];
    AdvanceParticleRestartOffset(local_mb_offsets[m + 1],
                                 static_cast<IOWrapperSizeT>(local_mb_counts[m]), 1);
  }
  const int local_npart =
      CheckedLocalParticleCount(local_mb_offsets[nmb_local], nrdata, nidata);
  ppart->nprtcl_thispack = local_npart;
  ValidateParticleRestartAllocation({nrdata, local_npart});
  Kokkos::realloc(ppart->prtcl_rdata, nrdata, local_npart);
  ValidateParticleRestartAllocation({nidata, local_npart});
  Kokkos::realloc(ppart->prtcl_idata, nidata, local_npart);

  std::vector<Real> packed_pr(static_cast<std::size_t>(local_npart)*nrdata, 0.0);
  std::vector<int> packed_pi(static_cast<std::size_t>(local_npart)*nidata, 0);

  const bool has_moments = (ref_has_moments > 0);
  const bool has_edge = (ref_has_edge > 0);

  HostArray5D<Real> h_mom("rst_mom", 1, 1, 1, 1, 1);
  if (has_moments) {
    if (ppart->moments.size() == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart expects particle moments, but moments are not allocated."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    h_mom = Kokkos::create_mirror_view(ppart->moments);
  }

  HostArray4D<Real> h_x1e("rst_jx1e", 1, 1, 1, 1);
  HostArray4D<Real> h_x2e("rst_jx2e", 1, 1, 1, 1);
  HostArray4D<Real> h_x3e("rst_jx3e", 1, 1, 1, 1);
  if (has_edge) {
    if (ppart->j_edge_x1e.size() == 0 || ppart->j_edge_x2e.size() == 0 ||
        ppart->j_edge_x3e.size() == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart expects particle edge-current state, but edge arrays are not "
                << "allocated." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    h_x1e = Kokkos::create_mirror_view(ppart->j_edge_x1e);
    h_x2e = Kokkos::create_mirror_view(ppart->j_edge_x2e);
    h_x3e = Kokkos::create_mirror_view(ppart->j_edge_x3e);
  }

  for (int r=0; r<meta.original_nranks; ++r) {
    auto &reqs = requests[r];
    if (reqs.empty() || !src_used[r]) continue;

    auto &sm = src_meta[r];
    IOWrapper srcfile;
    srcfile.Open(rank_paths[r].c_str(), IOWrapper::FileMode::read, true);

    for (const auto &req : reqs) {
      const int src_local = req.global_id - meta.gids_eachrank[r];
      if (src_local < 0 || src_local >= sm.nmb_section) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Restart metadata is inconsistent with source restart data in "
                  << "file '" << rank_paths[r] << "'." << std::endl;
        restart_utils::AbortOnFatalError();
      }

      const int cnt = local_mb_counts[req.local_index];
      const IOWrapperSizeT gstart = sm.mb_offsets[src_local];
      const IOWrapperSizeT lstart = local_mb_offsets[req.local_index];
      IOWrapperSizeT pr_off = sm.pr_real_offset;
      AdvanceParticleRestartOffset(pr_off, gstart,
                                   CheckedParticleRestartProduct(
                                       static_cast<IOWrapperSizeT>(nrdata),
                                       sizeof(Real)));
      IOWrapperSizeT pi_off = sm.pr_int_offset;
      AdvanceParticleRestartOffset(pi_off, gstart,
                                   CheckedParticleRestartProduct(
                                       static_cast<IOWrapperSizeT>(nidata), sizeof(int)));
      const IOWrapperSizeT pr_count =
          CheckedParticleRestartProduct(static_cast<IOWrapperSizeT>(cnt), nrdata);
      const IOWrapperSizeT pi_count =
          CheckedParticleRestartProduct(static_cast<IOWrapperSizeT>(cnt), nidata);
      const IOWrapperSizeT pr_start = CheckedParticleRestartProduct(lstart, nrdata);
      const IOWrapperSizeT pi_start = CheckedParticleRestartProduct(lstart, nidata);

      if (cnt > 0) {
        if (srcfile.Read_Reals_at(&(packed_pr[pr_start]), pr_count, pr_off,
                                  true) != static_cast<std::size_t>(pr_count)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to read particle restart real data from "
                    << "source restart file '" << rank_paths[r] << "'."
                    << std::endl;
          restart_utils::AbortOnFatalError();
        }
        if (srcfile.Read_bytes_at(&(packed_pi[pi_start]), sizeof(int), pi_count,
                                  pi_off, true) !=
            static_cast<std::size_t>(pi_count)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to read particle restart integer data from "
                    << "source restart file '" << rank_paths[r] << "'."
                    << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }

      if (has_moments) {
        auto mom_mb = Kokkos::subview(h_mom, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                      Kokkos::ALL, Kokkos::ALL);
        IOWrapperSizeT moff = sm.moments_offset;
        AdvanceParticleRestartOffset(moff, static_cast<IOWrapperSizeT>(src_local),
                                     static_cast<IOWrapperSizeT>(sm.moment_cnt)*
                                     sizeof(Real));
        if (srcfile.Read_Reals_at(mom_mb.data(), sm.moment_cnt, moff, true) !=
            static_cast<std::size_t>(sm.moment_cnt)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to read particle moment restart data from source "
                    << "restart file '" << rank_paths[r] << "'." << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }

      if (has_edge) {
        auto x1_mb = Kokkos::subview(h_x1e, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto x2_mb = Kokkos::subview(h_x2e, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto x3_mb = Kokkos::subview(h_x3e, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        IOWrapperSizeT x1off = sm.edge1_offset;
        AdvanceParticleRestartOffset(x1off, static_cast<IOWrapperSizeT>(src_local),
                                     static_cast<IOWrapperSizeT>(sm.edge1_cnt)*
                                     sizeof(Real));
        IOWrapperSizeT x2off = sm.edge2_offset;
        AdvanceParticleRestartOffset(x2off, static_cast<IOWrapperSizeT>(src_local),
                                     static_cast<IOWrapperSizeT>(sm.edge2_cnt)*
                                     sizeof(Real));
        IOWrapperSizeT x3off = sm.edge3_offset;
        AdvanceParticleRestartOffset(x3off, static_cast<IOWrapperSizeT>(src_local),
                                     static_cast<IOWrapperSizeT>(sm.edge3_cnt)*
                                     sizeof(Real));
        if (srcfile.Read_Reals_at(x1_mb.data(), sm.edge1_cnt, x1off, true) !=
              static_cast<std::size_t>(sm.edge1_cnt) ||
            srcfile.Read_Reals_at(x2_mb.data(), sm.edge2_cnt, x2off, true) !=
              static_cast<std::size_t>(sm.edge2_cnt) ||
            srcfile.Read_Reals_at(x3_mb.data(), sm.edge3_cnt, x3off, true) !=
              static_cast<std::size_t>(sm.edge3_cnt)) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl
                    << "Failed to read particle edge-current restart data from "
                    << "source restart file '" << rank_paths[r] << "'."
                    << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
    }

    srcfile.Close(true);
  }

  HostArray2D<Real> h_pr("rst_pr", nrdata, local_npart);
  HostArray2D<int> h_pi("rst_pi", nidata, local_npart);
  const int gids_local = pm->gids_eachrank[global_variable::my_rank];
  for (int p=0; p<local_npart; ++p) {
    for (int n=0; n<nrdata; ++n) {
      h_pr(n, p) = packed_pr[p*nrdata + n];
    }
    for (int n=0; n<nidata; ++n) {
      h_pi(n, p) = packed_pi[p*nidata + n];
    }
    ValidateRestoredParticleIDData(ppart, h_pi, p, gids_local, nmb_local);
  }
  Kokkos::deep_copy(ppart->prtcl_rdata, h_pr);
  Kokkos::deep_copy(ppart->prtcl_idata, h_pi);

  if (has_moments) {
    Kokkos::deep_copy(ppart->moments, h_mom);
  } else if (ppart->moments.size() > 0) {
    Kokkos::deep_copy(ppart->moments, static_cast<Real>(0.0));
  }

  if (has_edge) {
    Kokkos::deep_copy(ppart->j_edge_x1e, h_x1e);
    Kokkos::deep_copy(ppart->j_edge_x2e, h_x2e);
    Kokkos::deep_copy(ppart->j_edge_x3e, h_x3e);
  } else if (ppart->j_edge_x1e.size() > 0) {
    Kokkos::deep_copy(ppart->j_edge_x1e, static_cast<Real>(0.0));
    Kokkos::deep_copy(ppart->j_edge_x2e, static_cast<Real>(0.0));
    Kokkos::deep_copy(ppart->j_edge_x3e, static_cast<Real>(0.0));
  }

  pm->CountParticles();
}

void LoadSingleFileRestartData(Mesh *pm,
                               IOWrapperSizeT headeroffset,
                               IOWrapperSizeT data_stride,
                               int nout1, int nout2, int nout3,
                               int nhydro, int nmhd, int nrad,
                               int nforce, int nz4c, int nadm,
                               HostArray5D<Real> &ccin,
                               HostFaceFld4D<Real> &fcin) {
  MeshBlockPack *pack = pm->pmb_pack;
  int nmb = pack->nmb_thispack;
  hydro::Hydro* phydro = pack->phydro;
  mhd::MHD* pmhd = pack->pmhd;
  adm::ADM* padm = pack->padm;
  z4c::Z4c* pz4c = pack->pz4c;
  radiation::Radiation* prad = pack->prad;
  TurbulenceDriver* pturb = pack->pturb;

  const RestartMetaData &meta = pm->restart_meta;
  if (meta.file_name.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart metadata missing file name for single-file restart."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (meta.rank_eachmb.size() != static_cast<std::size_t>(pm->nmb_total)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Restart metadata inconsistent with MeshBlock count."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (meta.original_nranks <= 0 ||
      meta.gids_eachrank.size() != static_cast<std::size_t>(meta.original_nranks) ||
      meta.nmb_eachrank.size() != static_cast<std::size_t>(meta.original_nranks)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Restart metadata missing original rank layout."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  std::vector<std::vector<RestartBlockRequest>> requests(meta.original_nranks);
  for (int m=0; m<nmb; ++m) {
    int gid = pack->pmb->mb_gid.h_view(m);
    if (gid < 0 || gid >= pm->nmb_total) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Invalid MeshBlock gid encountered during restart."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    int src_rank = meta.rank_eachmb[gid];
    if (src_rank < 0 || src_rank >= meta.original_nranks) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Restart metadata contains invalid rank assignments."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    requests[src_rank].push_back({m, gid});
  }

  std::vector<std::string> rank_paths(meta.original_nranks);
  const std::string base_prefix = meta.base_dir.empty() ? "" :
      (meta.base_dir == "/" ? "/" : meta.base_dir + "/");
  for (int r=0; r<meta.original_nranks; ++r) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d", r);
    rank_paths[r] = base_prefix + rank_dir + "/" + meta.file_name;
  }

  const IOWrapperSizeT chunk_stride = data_stride;
  const int nout1p1 = CheckedParticleRestartExtentPlusOne(nout1);
  const int nout2p1 = CheckedParticleRestartExtentPlusOne(nout2);
  const int nout3p1 = CheckedParticleRestartExtentPlusOne(nout3);
  IOWrapperSizeT chunk_offset = 0;
  const IOWrapperSizeT hydro_offset = chunk_offset;
  AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1, nout2, nout3, nhydro});
  const IOWrapperSizeT mhd_cc_offset = chunk_offset;
  IOWrapperSizeT mhd_x1f_offset = chunk_offset;
  IOWrapperSizeT mhd_x2f_offset = chunk_offset;
  IOWrapperSizeT mhd_x3f_offset = chunk_offset;
  if (pmhd != nullptr) {
    AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1, nout2, nout3, nmhd});
    mhd_x1f_offset = chunk_offset;
    AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1p1, nout2, nout3});
    mhd_x2f_offset = chunk_offset;
    AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1, nout2p1, nout3});
    mhd_x3f_offset = chunk_offset;
    AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1, nout2, nout3p1});
  }
  const IOWrapperSizeT rad_offset = chunk_offset;
  AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1, nout2, nout3, nrad});
  const IOWrapperSizeT turb_offset = chunk_offset;
  if (pturb != nullptr && nforce > 0) {
    AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1, nout2, nout3, nforce});
  }
  const IOWrapperSizeT z4c_adm_offset = chunk_offset;
  if (pz4c != nullptr) {
    AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1, nout2, nout3, nz4c});
  } else if (padm != nullptr) {
    AdvanceParticleRestartRealArrayOffset(chunk_offset, {nout1, nout2, nout3, nadm});
  }
  if (chunk_offset != chunk_stride) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Restart data chunk size mismatch, restart file is broken."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  auto chunk_base = [&](int src_rank, int global_id) -> IOWrapperSizeT {
    int start_gid = meta.gids_eachrank[src_rank];
    int local_index = global_id - start_gid;
    if (local_index < 0 || local_index >= meta.nmb_eachrank[src_rank]) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Restart metadata inconsistent with MeshBlock ids."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    IOWrapperSizeT offset = headeroffset;
    AdvanceParticleRestartOffset(offset, static_cast<IOWrapperSizeT>(local_index),
                                 chunk_stride);
    return offset;
  };

  if (phydro != nullptr && nhydro > 0) {
    ValidateParticleRestartAllocation({nmb, nhydro, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nhydro, nout3, nout2, nout1);
    for (int r=0; r<meta.original_nranks; ++r) {
      auto &reqs = requests[r];
      if (reqs.empty()) continue;
      IOWrapper srcfile;
      srcfile.Open(rank_paths[r].c_str(), IOWrapper::FileMode::read, true);
      for (const auto &req : reqs) {
        auto mbptr = Kokkos::subview(ccin, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (mbcnt > 0) {
          IOWrapperSizeT base = chunk_base(r, req.global_id);
          if (srcfile.Read_Reals_at(
                  mbptr.data(), mbcnt,
                  CheckedParticleRestartDataOffset(base, hydro_offset), true)
              != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "CC hydro data not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }
      }
      srcfile.Close(true);
    }
    Kokkos::deep_copy(Kokkos::subview(phydro->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
  }

  if (pmhd != nullptr && nmhd > 0) {
    ValidateParticleRestartAllocation({nmb, nmhd, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nmhd, nout3, nout2, nout1);
    ValidateParticleRestartAllocation({nmb, nout3, nout2, nout1p1});
    Kokkos::realloc(fcin.x1f, nmb, nout3, nout2, nout1p1);
    ValidateParticleRestartAllocation({nmb, nout3, nout2p1, nout1});
    Kokkos::realloc(fcin.x2f, nmb, nout3, nout2p1, nout1);
    ValidateParticleRestartAllocation({nmb, nout3p1, nout2, nout1});
    Kokkos::realloc(fcin.x3f, nmb, nout3p1, nout2, nout1);
    for (int r=0; r<meta.original_nranks; ++r) {
      auto &reqs = requests[r];
      if (reqs.empty()) continue;
      IOWrapper srcfile;
      srcfile.Open(rank_paths[r].c_str(), IOWrapper::FileMode::read, true);
      for (const auto &req : reqs) {
        IOWrapperSizeT base = chunk_base(r, req.global_id);
        auto mbptr = Kokkos::subview(ccin, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (mbcnt > 0) {
          if (srcfile.Read_Reals_at(
                  mbptr.data(), mbcnt,
                  CheckedParticleRestartDataOffset(base, mhd_cc_offset), true)
              != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "CC mhd data not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }

        auto x1fptr = Kokkos::subview(fcin.x1f, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                       Kokkos::ALL);
        auto fldcnt = x1fptr.size();
        if (fldcnt > 0) {
          if (srcfile.Read_Reals_at(
                  x1fptr.data(), fldcnt,
                  CheckedParticleRestartDataOffset(base, mhd_x1f_offset), true)
              != fldcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl
                      << "Input b0.x1f field not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }

        auto x2fptr = Kokkos::subview(fcin.x2f, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                       Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (fldcnt > 0) {
          if (srcfile.Read_Reals_at(
                  x2fptr.data(), fldcnt,
                  CheckedParticleRestartDataOffset(base, mhd_x2f_offset), true)
              != fldcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl
                      << "Input b0.x2f field not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }

        auto x3fptr = Kokkos::subview(fcin.x3f, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                       Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (fldcnt > 0) {
          if (srcfile.Read_Reals_at(
                  x3fptr.data(), fldcnt,
                  CheckedParticleRestartDataOffset(base, mhd_x3f_offset), true)
              != fldcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl
                      << "Input b0.x3f field not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }
      }
      srcfile.Close(true);
    }
    Kokkos::deep_copy(Kokkos::subview(pmhd->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x1f);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x2f);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x3f);
  }

  if (prad != nullptr && nrad > 0) {
    ValidateParticleRestartAllocation({nmb, nrad, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nrad, nout3, nout2, nout1);
    for (int r=0; r<meta.original_nranks; ++r) {
      auto &reqs = requests[r];
      if (reqs.empty()) continue;
      IOWrapper srcfile;
      srcfile.Open(rank_paths[r].c_str(), IOWrapper::FileMode::read, true);
      for (const auto &req : reqs) {
        auto mbptr = Kokkos::subview(ccin, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (mbcnt > 0) {
          IOWrapperSizeT base = chunk_base(r, req.global_id);
          if (srcfile.Read_Reals_at(
                  mbptr.data(), mbcnt,
                  CheckedParticleRestartDataOffset(base, rad_offset), true)
              != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "CC rad data not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }
      }
      srcfile.Close(true);
    }
    Kokkos::deep_copy(Kokkos::subview(prad->i0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
  }

  if (pturb != nullptr && nforce > 0) {
    ValidateParticleRestartAllocation({nmb, nforce, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nforce, nout3, nout2, nout1);
    for (int r=0; r<meta.original_nranks; ++r) {
      auto &reqs = requests[r];
      if (reqs.empty()) continue;
      IOWrapper srcfile;
      srcfile.Open(rank_paths[r].c_str(), IOWrapper::FileMode::read, true);
      for (const auto &req : reqs) {
        auto mbptr = Kokkos::subview(ccin, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (mbcnt > 0) {
          IOWrapperSizeT base = chunk_base(r, req.global_id);
          if (srcfile.Read_Reals_at(
                  mbptr.data(), mbcnt,
                  CheckedParticleRestartDataOffset(base, turb_offset), true)
              != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "CC turb data not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }
      }
      srcfile.Close(true);
    }
    Kokkos::deep_copy(Kokkos::subview(pturb->force, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
  }

  if (pz4c != nullptr && nz4c > 0) {
    ValidateParticleRestartAllocation({nmb, nz4c, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nz4c, nout3, nout2, nout1);
    for (int r=0; r<meta.original_nranks; ++r) {
      auto &reqs = requests[r];
      if (reqs.empty()) continue;
      IOWrapper srcfile;
      srcfile.Open(rank_paths[r].c_str(), IOWrapper::FileMode::read, true);
      for (const auto &req : reqs) {
        auto mbptr = Kokkos::subview(ccin, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (mbcnt > 0) {
          IOWrapperSizeT base = chunk_base(r, req.global_id);
          if (srcfile.Read_Reals_at(
                  mbptr.data(), mbcnt,
                  CheckedParticleRestartDataOffset(base, z4c_adm_offset), true)
              != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "CC z4c data not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }
      }
      srcfile.Close(true);
    }
    Kokkos::deep_copy(Kokkos::subview(pz4c->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    pz4c->Z4cToADM(pm->pmb_pack);
  } else if (padm != nullptr && nadm > 0) {
    ValidateParticleRestartAllocation({nmb, nadm, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nadm, nout3, nout2, nout1);
    for (int r=0; r<meta.original_nranks; ++r) {
      auto &reqs = requests[r];
      if (reqs.empty()) continue;
      IOWrapper srcfile;
      srcfile.Open(rank_paths[r].c_str(), IOWrapper::FileMode::read, true);
      for (const auto &req : reqs) {
        auto mbptr = Kokkos::subview(ccin, req.local_index, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (mbcnt > 0) {
          IOWrapperSizeT base = chunk_base(r, req.global_id);
          if (srcfile.Read_Reals_at(
                  mbptr.data(), mbcnt,
                  CheckedParticleRestartDataOffset(base, z4c_adm_offset), true)
              != mbcnt) {
            std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                      << std::endl << "CC adm data not read correctly from rst file, "
                      << "restart file is broken." << std::endl;
            restart_utils::AbortOnFatalError();
          }
        }
      }
      srcfile.Close(true);
    }
    Kokkos::deep_copy(Kokkos::subview(padm->u_adm, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
  }
}

void LoadParticleRestartData(Mesh *pm,
                             IOWrapper &resfile,
                             bool single_file_per_rank,
                             IOWrapperSizeT headeroffset,
                             IOWrapperSizeT data_stride,
                             int nout1, int nout2, int nout3) {
  auto *ppart = pm->pmb_pack->ppart;
  if (ppart == nullptr) return;

  if (single_file_per_rank) {
    LoadParticleRestartDataSingleFile(pm, headeroffset, data_stride,
                                      nout1, nout2, nout3);
    return;
  }

  constexpr std::uint64_t kPicMagic = 0x5049435253543031ULL;
  constexpr int kPicVersion = particles::Particles::PIC_RESTART_SCHEMA_VERSION;

  IOWrapperSizeT section_offset = headeroffset;
  AdvanceParticleRestartOffset(section_offset, static_cast<IOWrapperSizeT>(pm->nmb_total),
                               data_stride);

  std::uint64_t pic_magic = 0;
  int has_pic_section = 0;
  if (global_variable::my_rank == 0) {
    if (resfile.Read_bytes_at(&pic_magic, 1, sizeof(std::uint64_t), section_offset,
                              false) == sizeof(std::uint64_t)) {
      has_pic_section = 1;
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(&has_pic_section, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&pic_magic, sizeof(std::uint64_t), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  if (!has_pic_section) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart state is missing from restart file."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (pic_magic != kPicMagic) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Restart file has an unknown particle restart marker."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  IOWrapperSizeT rd_offset = section_offset;
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(std::uint64_t));
  int version = -1;
  int nmb_section = 0;
  int nrdata = 0;
  int nidata = 0;
  int rst_nout1 = 0;
  int rst_nout2 = 0;
  int rst_nout3 = 0;
  int has_moments = 0;
  int has_edge = 0;
  int moment_cnt = 0;
  int edge1_cnt = 0;
  int edge2_cnt = 0;
  int edge3_cnt = 0;
  int state_kind = -1;
  int physical_mode = -1;
  Real cr_light_speed = 0.0;
  std::array<int, particles::Particles::NPIC_RESTART_MODEL_INTS> model_ints;
  std::array<Real, particles::Particles::NPIC_RESTART_MODEL_REALS> model_reals;
  IOWrapperSizeT npart_section = 0;

  auto read_int_meta = [&](int &val, IOWrapperSizeT off, const char *name) {
    if (global_variable::my_rank == 0) {
      if (resfile.Read_bytes_at(&val, sizeof(int), 1, off, false) != 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to read particle restart metadata field '" << name << "'."
                  << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(&val, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  };

  read_int_meta(version, rd_offset, "version");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  if (version != kPicVersion) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Unsupported particle restart version " << version << "."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  read_int_meta(nmb_section, rd_offset, "nmb_section");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(nrdata, rd_offset, "nrdata");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(nidata, rd_offset, "nidata");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(rst_nout1, rd_offset, "nout1");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(rst_nout2, rd_offset, "nout2");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(rst_nout3, rd_offset, "nout3");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(has_moments, rd_offset, "has_moments");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(has_edge, rd_offset, "has_edge");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(moment_cnt, rd_offset, "moment_cnt");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(edge1_cnt, rd_offset, "edge1_cnt");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(edge2_cnt, rd_offset, "edge2_cnt");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(edge3_cnt, rd_offset, "edge3_cnt");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(state_kind, rd_offset, "state_kind");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  read_int_meta(physical_mode, rd_offset, "physical_mode");
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(int));
  if (global_variable::my_rank == 0) {
    if (resfile.Read_Reals_at(&cr_light_speed, 1, rd_offset, false) != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart metadata field 'cr_light_speed'."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(&cr_light_speed, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
#endif
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(Real));
  if (global_variable::my_rank == 0) {
    if (resfile.Read_bytes_at(model_ints.data(), sizeof(int), model_ints.size(),
                              rd_offset, false) != model_ints.size()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart extension-model integer metadata."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(model_ints.data(), model_ints.size(), MPI_INT, 0, MPI_COMM_WORLD);
#endif
  AdvanceParticleRestartOffset(rd_offset, model_ints.size(), sizeof(int));
  if (global_variable::my_rank == 0) {
    if (resfile.Read_Reals_at(model_reals.data(), model_reals.size(), rd_offset, false)
          != model_reals.size()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart extension-model numerical metadata."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(model_reals.data(), model_reals.size(), MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
#endif
  AdvanceParticleRestartOffset(rd_offset, model_reals.size(), sizeof(Real));

  if (global_variable::my_rank == 0) {
    if (resfile.Read_bytes_at(&npart_section, sizeof(IOWrapperSizeT), 1, rd_offset,
                              false) != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart count metadata."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(&npart_section, sizeof(IOWrapperSizeT), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  AdvanceParticleRestartOffset(rd_offset, 1, sizeof(IOWrapperSizeT));

  if (nmb_section != pm->nmb_total) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart MeshBlock count mismatch (file=" << nmb_section
              << ", runtime=" << pm->nmb_total << ")." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (nrdata <= 0 || nidata <= 0 ||
      nrdata != ppart->nrdata || nidata != ppart->nidata) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart data layout mismatch." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const int expected_state_kind = ppart->UsesRelativisticCRState() ? 1 : 0;
  const int expected_physical_mode = static_cast<int>(ppart->pic_physical_mode);
  if (state_kind != expected_state_kind ||
      physical_mode != expected_physical_mode ||
      cr_light_speed != ppart->pic_cr_light_speed ||
      !ppart->MatchesRestartModelMetadata(model_ints, model_reals) ||
      !ppart->RestoreRestartModelState(model_reals)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart physical-model metadata mismatch." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (rst_nout1 != nout1 || rst_nout2 != nout2 || rst_nout3 != nout3) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart mesh extents mismatch for moment arrays."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if ((has_moments != 0 && has_moments != 1) || (has_edge != 0 && has_edge != 1)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart optional-array flags are invalid."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  const int expected_moment_cnt = CheckedParticleRestartIntProduct(
      {particles::Particles::NMOM, nout3, nout2, nout1});
  const int nout1p1 = CheckedParticleRestartExtentPlusOne(nout1);
  const int nout2p1 = CheckedParticleRestartExtentPlusOne(nout2);
  const int nout3p1 = CheckedParticleRestartExtentPlusOne(nout3);
  const int expected_edge1_cnt =
      CheckedParticleRestartIntProduct({nout3p1, nout2p1, nout1});
  const int expected_edge2_cnt =
      CheckedParticleRestartIntProduct({nout3p1, nout2, nout1p1});
  const int expected_edge3_cnt =
      CheckedParticleRestartIntProduct({nout3, nout2p1, nout1p1});
  if (moment_cnt != expected_moment_cnt) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart moment size mismatch." << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (has_edge != 0 &&
      (edge1_cnt != expected_edge1_cnt || edge2_cnt != expected_edge2_cnt ||
       edge3_cnt != expected_edge3_cnt)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart edge-current size mismatch." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  const IOWrapperSizeT mb_count_offset = rd_offset;
  std::vector<int> mb_counts(nmb_section, 0);
  if (global_variable::my_rank == 0 && nmb_section > 0) {
    if (resfile.Read_bytes_at(mb_counts.data(), sizeof(int), nmb_section, mb_count_offset,
                              false) != static_cast<std::size_t>(nmb_section)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart MeshBlock counts."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  if (nmb_section > 0) {
    MPI_Bcast(mb_counts.data(), nmb_section, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif

  std::vector<IOWrapperSizeT> mb_offsets(nmb_section + 1, 0);
  for (int m=0; m<nmb_section; ++m) {
    if (mb_counts[m] < 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Particle restart MeshBlock count table contains a negative count."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    mb_offsets[m + 1] = mb_offsets[m];
    AdvanceParticleRestartOffset(mb_offsets[m + 1],
                                 static_cast<IOWrapperSizeT>(mb_counts[m]), 1);
  }
  if (mb_offsets[nmb_section] != npart_section) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart count table is inconsistent." << std::endl;
    restart_utils::AbortOnFatalError();
  }

  IOWrapperSizeT data_offset = mb_count_offset;
  AdvanceParticleRestartOffset(data_offset, static_cast<IOWrapperSizeT>(nmb_section),
                               sizeof(int));
  const IOWrapperSizeT pr_real_offset = data_offset;
  AdvanceParticleRestartOffset(data_offset, npart_section,
                               CheckedParticleRestartProduct(
                                   static_cast<IOWrapperSizeT>(nrdata), sizeof(Real)));
  const IOWrapperSizeT pr_int_offset = data_offset;
  AdvanceParticleRestartOffset(data_offset, npart_section,
                               CheckedParticleRestartProduct(
                                   static_cast<IOWrapperSizeT>(nidata), sizeof(int)));
  const IOWrapperSizeT moments_offset = data_offset;
  if (has_moments) {
    AdvanceParticleRestartOffset(data_offset, static_cast<IOWrapperSizeT>(nmb_section),
                                 CheckedParticleRestartProduct(
                                     static_cast<IOWrapperSizeT>(moment_cnt),
                                     sizeof(Real)));
  }
  const IOWrapperSizeT edge1_offset = data_offset;
  if (has_edge) {
    AdvanceParticleRestartOffset(data_offset, static_cast<IOWrapperSizeT>(nmb_section),
                                 CheckedParticleRestartProduct(
                                     static_cast<IOWrapperSizeT>(edge1_cnt),
                                     sizeof(Real)));
  }
  const IOWrapperSizeT edge2_offset = data_offset;
  if (has_edge) {
    AdvanceParticleRestartOffset(data_offset, static_cast<IOWrapperSizeT>(nmb_section),
                                 CheckedParticleRestartProduct(
                                     static_cast<IOWrapperSizeT>(edge2_cnt),
                                     sizeof(Real)));
  }
  const IOWrapperSizeT edge3_offset = data_offset;
  if (has_edge) {
    AdvanceParticleRestartOffset(data_offset, static_cast<IOWrapperSizeT>(nmb_section),
                                 CheckedParticleRestartProduct(
                                     static_cast<IOWrapperSizeT>(edge3_cnt),
                                     sizeof(Real)));
  }
  if (data_offset > resfile.GetSize(false)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "Particle restart section exceeds source artifact bounds."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }

  const int nmb_local = pm->nmb_thisrank;
  const int gids_local = pm->gids_eachrank[global_variable::my_rank];
  std::vector<int> local_mb_counts(nmb_local, 0);
  for (int m=0; m<nmb_local; ++m) {
    local_mb_counts[m] = mb_counts[gids_local + m];
  }
  std::vector<IOWrapperSizeT> local_mb_offsets(nmb_local + 1, 0);
  for (int m=0; m<nmb_local; ++m) {
    local_mb_offsets[m + 1] = local_mb_offsets[m];
    AdvanceParticleRestartOffset(local_mb_offsets[m + 1],
                                 static_cast<IOWrapperSizeT>(local_mb_counts[m]), 1);
  }

  const int local_npart =
      CheckedLocalParticleCount(local_mb_offsets[nmb_local], nrdata, nidata);
  ppart->nprtcl_thispack = local_npart;
  ValidateParticleRestartAllocation({nrdata, local_npart});
  Kokkos::realloc(ppart->prtcl_rdata, nrdata, local_npart);
  ValidateParticleRestartAllocation({nidata, local_npart});
  Kokkos::realloc(ppart->prtcl_idata, nidata, local_npart);

  std::vector<Real> packed_pr(static_cast<std::size_t>(local_npart)*nrdata, 0.0);
  std::vector<int> packed_pi(static_cast<std::size_t>(local_npart)*nidata, 0);

  for (int m=0; m<nmb_local; ++m) {
    const int cnt = local_mb_counts[m];
    if (cnt <= 0) continue;
    const int gid = gids_local + m;
    const IOWrapperSizeT gstart = mb_offsets[gid];
    const IOWrapperSizeT lstart = local_mb_offsets[m];
    IOWrapperSizeT pr_off = pr_real_offset;
    AdvanceParticleRestartOffset(pr_off, gstart,
                                 CheckedParticleRestartProduct(
                                     static_cast<IOWrapperSizeT>(nrdata), sizeof(Real)));
    IOWrapperSizeT pi_off = pr_int_offset;
    AdvanceParticleRestartOffset(pi_off, gstart,
                                 CheckedParticleRestartProduct(
                                     static_cast<IOWrapperSizeT>(nidata), sizeof(int)));
    const IOWrapperSizeT pr_count =
        CheckedParticleRestartProduct(static_cast<IOWrapperSizeT>(cnt), nrdata);
    const IOWrapperSizeT pi_count =
        CheckedParticleRestartProduct(static_cast<IOWrapperSizeT>(cnt), nidata);
    const IOWrapperSizeT pr_start = CheckedParticleRestartProduct(lstart, nrdata);
    const IOWrapperSizeT pi_start = CheckedParticleRestartProduct(lstart, nidata);
    if (resfile.Read_Reals_at(&(packed_pr[pr_start]), pr_count, pr_off, false) !=
        static_cast<std::size_t>(pr_count)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart real data." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    if (resfile.Read_bytes_at(&(packed_pi[pi_start]), sizeof(int), pi_count,
                              pi_off, false) !=
        static_cast<std::size_t>(pi_count)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Failed to read particle restart integer data." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }

  HostArray2D<Real> h_pr("rst_pr", nrdata, local_npart);
  HostArray2D<int> h_pi("rst_pi", nidata, local_npart);
  for (int p=0; p<local_npart; ++p) {
    for (int n=0; n<nrdata; ++n) {
      h_pr(n, p) = packed_pr[p*nrdata + n];
    }
    for (int n=0; n<nidata; ++n) {
      h_pi(n, p) = packed_pi[p*nidata + n];
    }
    ValidateRestoredParticleIDData(ppart, h_pi, p, gids_local, nmb_local);
  }
  Kokkos::deep_copy(ppart->prtcl_rdata, h_pr);
  Kokkos::deep_copy(ppart->prtcl_idata, h_pi);

  if (has_moments) {
    if (ppart->moments.size() == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart expects particle moments, but moments are not allocated."
                << std::endl;
      restart_utils::AbortOnFatalError();
    }
    auto h_mom = Kokkos::create_mirror_view(ppart->moments);
    for (int m=0; m<nmb_local; ++m) {
      const int gid = gids_local + m;
      auto mom_mb = Kokkos::subview(h_mom, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                    Kokkos::ALL);
      IOWrapperSizeT moff = moments_offset;
      AdvanceParticleRestartOffset(moff, static_cast<IOWrapperSizeT>(gid),
                                   CheckedParticleRestartProduct(
                                       static_cast<IOWrapperSizeT>(moment_cnt),
                                       sizeof(Real)));
      if (resfile.Read_Reals_at(mom_mb.data(), moment_cnt, moff, false) !=
          static_cast<std::size_t>(moment_cnt)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to read particle moment restart data." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
    Kokkos::deep_copy(ppart->moments, h_mom);
  } else if (ppart->moments.size() > 0) {
    Kokkos::deep_copy(ppart->moments, static_cast<Real>(0.0));
  }

  if (has_edge) {
    if (ppart->j_edge_x1e.size() == 0 || ppart->j_edge_x2e.size() == 0 ||
        ppart->j_edge_x3e.size() == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "Restart expects particle edge-current state, but edge arrays are not "
                << "allocated." << std::endl;
      restart_utils::AbortOnFatalError();
    }
    auto h_x1e = Kokkos::create_mirror_view(ppart->j_edge_x1e);
    auto h_x2e = Kokkos::create_mirror_view(ppart->j_edge_x2e);
    auto h_x3e = Kokkos::create_mirror_view(ppart->j_edge_x3e);
    for (int m=0; m<nmb_local; ++m) {
      const int gid = gids_local + m;
      auto x1_mb = Kokkos::subview(h_x1e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      auto x2_mb = Kokkos::subview(h_x2e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      auto x3_mb = Kokkos::subview(h_x3e, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
      IOWrapperSizeT x1off = edge1_offset;
      AdvanceParticleRestartOffset(x1off, static_cast<IOWrapperSizeT>(gid),
                                   CheckedParticleRestartProduct(
                                       static_cast<IOWrapperSizeT>(edge1_cnt),
                                       sizeof(Real)));
      IOWrapperSizeT x2off = edge2_offset;
      AdvanceParticleRestartOffset(x2off, static_cast<IOWrapperSizeT>(gid),
                                   CheckedParticleRestartProduct(
                                       static_cast<IOWrapperSizeT>(edge2_cnt),
                                       sizeof(Real)));
      IOWrapperSizeT x3off = edge3_offset;
      AdvanceParticleRestartOffset(x3off, static_cast<IOWrapperSizeT>(gid),
                                   CheckedParticleRestartProduct(
                                       static_cast<IOWrapperSizeT>(edge3_cnt),
                                       sizeof(Real)));
      if (resfile.Read_Reals_at(x1_mb.data(), edge1_cnt, x1off, false) !=
            static_cast<std::size_t>(edge1_cnt) ||
          resfile.Read_Reals_at(x2_mb.data(), edge2_cnt, x2off, false) !=
            static_cast<std::size_t>(edge2_cnt) ||
          resfile.Read_Reals_at(x3_mb.data(), edge3_cnt, x3off, false) !=
            static_cast<std::size_t>(edge3_cnt)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "Failed to read particle edge-current restart data." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
    Kokkos::deep_copy(ppart->j_edge_x1e, h_x1e);
    Kokkos::deep_copy(ppart->j_edge_x2e, h_x2e);
    Kokkos::deep_copy(ppart->j_edge_x3e, h_x3e);
  } else {
    if (ppart->j_edge_x1e.size() > 0) {
      Kokkos::deep_copy(ppart->j_edge_x1e, static_cast<Real>(0.0));
      Kokkos::deep_copy(ppart->j_edge_x2e, static_cast<Real>(0.0));
      Kokkos::deep_copy(ppart->j_edge_x3e, static_cast<Real>(0.0));
    }
  }

  pm->CountParticles();
}

}  // namespace

//----------------------------------------------------------------------------------------
// default constructor, calls pgen function.

ProblemGenerator::ProblemGenerator(ParameterInput *pin, Mesh *pm) :
    user_bcs(false),
    user_srcs(false),
    user_hist(false),
    user_work_in_loop(false),
    pmy_mesh_(pm) {
  // check for user-defined boundary conditions
  for (int dir=0; dir<6; ++dir) {
    if (pm->mesh_bcs[dir] == BoundaryFlag::user) {
      user_bcs = true;
    }
  }

  user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  user_hist = pin->GetOrAddBoolean("problem","user_hist",false);
  user_work_in_loop = pin->GetOrAddBoolean("problem", "user_work_in_loop", false);

#if USER_PROBLEM_ENABLED
  // call user-defined problem generator
  UserProblem(pin, false);
#else
  // else read name of built-in pgen from <problem> block in input file, and call
  std::string pgen_fun_name = pin->GetOrAddString("problem", "pgen_name", "none");

  if (pgen_fun_name.compare("advection") == 0) {
    Advection(pin, false);
  } else if (pgen_fun_name.compare("cpaw") == 0) {
    AlfvenWave(pin, false);
  } else if (pgen_fun_name.compare("gr_bondi") == 0) {
    BondiAccretion(pin, false);
  } else if (pgen_fun_name.compare("tetrad") == 0) {
    CheckOrthonormalTetrad(pin, false);
  } else if (pgen_fun_name.compare("hohlraum") == 0) {
    Hohlraum(pin, false);
  } else if (pgen_fun_name.compare("linear_wave") == 0) {
    LinearWave(pin, false);
  } else if (pgen_fun_name.compare("implode") == 0) {
    LWImplode(pin, false);
  } else if (pgen_fun_name.compare("gr_monopole") == 0) {
    Monopole(pin, false);
  } else if (pgen_fun_name.compare("orszag_tang") == 0) {
    OrszagTang(pin, false);
  } else if (pgen_fun_name.compare("pic_paper_smooth_tsc_interface") == 0) {
    PICPaperSmoothTSCInterface(pin, false);
  } else if (pgen_fun_name.compare("pic_parallel_shock") == 0) {
    PICParallelShock(pin, false);
  } else if (pgen_fun_name.compare("q006_paper_multispecies_oscillation") == 0) {
    Q006PaperMultispeciesOscillation(pin, false);
  } else if (
      pgen_fun_name.compare("q006_paper_multispecies_oscillation_runtime_local")
      == 0) {
    Q006PaperMultispeciesOscillationRuntimeLocal(pin, false);
  } else if (pgen_fun_name.compare("q007_paper_crsi_linear_preparation") == 0) {
    Q007PaperCRSILinearPreparation(pin, false);
  } else if (pgen_fun_name.compare("q007_paper_crpai_linear_preparation") == 0) {
    Q007PaperCRPAILinearPreparation(pin, false);
  } else if (pgen_fun_name.compare("q023_paper_bell_linear") == 0) {
    Q023PaperBellLinear(pin, false);
  } else if (pgen_fun_name.compare("q043_bell_current_volume_aware") == 0) {
    Q043BellCurrentVolumeAware(pin, false);
  } else if (pgen_fun_name.compare("q029_hall_bell_linear") == 0) {
    Q029HallBellLinear(pin, false);
  } else if (pgen_fun_name.compare("q032_reduced_static_neutral_local") == 0) {
    Q032ReducedStaticNeutralLocal(pin, false);
  } else if (pgen_fun_name.compare("q033_crpai_transport_runtime_local") == 0) {
    Q033CRPAITransportRuntimeLocal(pin, false);
  } else if (pgen_fun_name.compare("rad_linear_wave") == 0) {
    RadiationLinearWave(pin, false);
  } else if (pgen_fun_name.compare("shock_tube") == 0) {
    ShockTube(pin, false);
  } else if (pgen_fun_name.compare("z4c_linear_wave") == 0) {
    Z4cLinearWave(pin, false);
  } else if (pgen_fun_name.compare("spherical_collapse") == 0) {
    SphericalCollapse(pin, false);
  } else if (pgen_fun_name.compare("diffusion") == 0) {
    Diffusion(pin, false);
  // else, name not set on command line or input file, print warning and quit
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "Problem generator name could not be found in <problem> block in input file"
        << std::endl
        << "and it was not set by -D PROBLEM option on cmake command line during build"
        << std::endl
        << "Rerun cmake with -D PROBLEM=file to specify custom problem generator file"
        << std::endl;;
    restart_utils::AbortOnFatalError();
  }
#endif

  // Check that user defined BCs were enrolled if needed
  if (user_bcs) {
    if (user_bcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User BCs specified in <mesh> block, but not enrolled "
                << "by SetProblemData()." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  // Check that user defined srcterms were enrolled if needed
  if (user_srcs) {
    if (user_srcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User SRCs specified in <problem> block, but not "
                << "enrolled by UserProblem()." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  // Check that user defined history outputs were enrolled if needed
  if (user_hist) {
    if (user_hist_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User history output specified in <problem> block, but "
                << "not enrolled by UserProblem()." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  if (user_work_in_loop) {
    if (user_work_in_loop_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User work-in-loop callback specified in <problem> "
                << "block, but not enrolled by UserProblem()." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
}

//----------------------------------------------------------------------------------------
// constructor for restarts
// When called, data needed to rebuild mesh has been read from restart file by
// Mesh::BuildTreeFromRestart() function. This constructor reads from the restart file and
// initializes all the dependent variables (u0,b0,etc) stored in each Physics class. It
// also calls ProblemGenerator::SetProblemData() function to set any user-defined BCs,
// and any data necessary for restart runs to continue correctly.

ProblemGenerator::ProblemGenerator(ParameterInput *pin, Mesh *pm, IOWrapper resfile,
                                   bool single_file_per_rank) :
    user_bcs(false),
    user_srcs(false),
    user_hist(false),
    user_work_in_loop(false),
    pmy_mesh_(pm) {
  // check for user-defined boundary conditions
  for (int dir=0; dir<6; ++dir) {
    if (pm->mesh_bcs[dir] == BoundaryFlag::user) {
      user_bcs = true;
    }
  }
  user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  user_hist = pin->GetOrAddBoolean("problem","user_hist",false);
  user_work_in_loop = pin->GetOrAddBoolean("problem", "user_work_in_loop", false);

  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = CheckedParticleRestartOutputExtent(indcs.nx1, indcs.ng);
  int nout2 = CheckedParticleRestartOutputExtent(indcs.nx2, indcs.ng);
  int nout3 = CheckedParticleRestartOutputExtent(indcs.nx3, indcs.ng);
  const int nout1p1 = CheckedParticleRestartExtentPlusOne(nout1);
  const int nout2p1 = CheckedParticleRestartExtentPlusOne(nout2);
  const int nout3p1 = CheckedParticleRestartExtentPlusOne(nout3);
  int nmb = pm->pmb_pack->nmb_thispack;
  // calculate total number of CC variables
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  adm::ADM* padm = pm->pmb_pack->padm;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  radiation::Radiation* prad=pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  int nrad = 0, nhydro = 0, nmhd = 0, nforce = 3, nadm = 0, nz4c = 0;
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }

  // root process reads z4c last_output_time and tracker data
  if (pz4c != nullptr) {
    Real last_output_time;
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_Reals(&last_output_time, 1,single_file_per_rank) != 1) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "z4c::last_output_time data size read from restart "
                  << "file is incorrect, restart file is broken." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&last_output_time, sizeof(Real), MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    pz4c->last_output_time = last_output_time;

    for (auto &pt : pz4c->ptracker) {
      Real pos[3];
      if (global_variable::my_rank == 0 || single_file_per_rank) {
        if (resfile.Read_Reals(&pos[0], 3, single_file_per_rank) != 3) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "compact object tracker data size read from restart "
                    << "file is incorrect, restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
      }
#if MPI_PARALLEL_ENABLED
      if (!single_file_per_rank) {
        MPI_Bcast(&pos[0], 3*sizeof(Real), MPI_CHAR, 0, MPI_COMM_WORLD);
      }
#endif
      pt.SetPos(&pos[0]);
    }
  }

  if (pturb != nullptr) {
    // root process reads size the random seed
    char rng_data[sizeof(RNG_State)];
    // the master process reads the variables data
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      if (resfile.Read_bytes(rng_data, 1, sizeof(RNG_State), single_file_per_rank)
          != sizeof(RNG_State)) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "RNG data size read from restart file is incorrect, "
                  << "restart file is broken." << std::endl;
        restart_utils::AbortOnFatalError();
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      // then broadcast the RNG information
      MPI_Bcast(rng_data, sizeof(RNG_State), MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
    std::memcpy(&(pturb->rstate), &(rng_data[0]), sizeof(RNG_State));
  }

  // root process reads size of CC and FC data arrays from restart file
  IOWrapperSizeT variablesize = sizeof(IOWrapperSizeT);
  char variabledata[sizeof(IOWrapperSizeT)];
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read_bytes(variabledata, 1, variablesize, single_file_per_rank)
        != variablesize) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Variable data size read from restart file is incorrect, "
                << "restart file is broken." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
#if MPI_PARALLEL_ENABLED
  // then broadcast the datasize information
  if (!single_file_per_rank) {
    MPI_Bcast(variabledata, variablesize, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif
  IOWrapperSizeT data_size;
  std::memcpy(&data_size, &(variabledata[0]), sizeof(IOWrapperSizeT));

  // calculate total number of CC variables
  IOWrapperSizeT headeroffset;
  // master process gets file offset
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    headeroffset = resfile.GetPosition(single_file_per_rank);
  }
#if MPI_PARALLEL_ENABLED
  // then broadcasts it
  if (!single_file_per_rank) {
    MPI_Bcast(&headeroffset, sizeof(IOWrapperSizeT), MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif

  IOWrapperSizeT data_size_ = 0;
  if (phydro != nullptr) {
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1, nout2, nout3, nhydro});
  }
  if (pmhd != nullptr) {
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1, nout2, nout3, nmhd});
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1p1, nout2, nout3});
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1, nout2p1, nout3});
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1, nout2, nout3p1});
  }
  if (prad != nullptr) {
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1, nout2, nout3, nrad});
  }
  if (pturb != nullptr) {
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1, nout2, nout3, nforce});
  }
  if (pz4c != nullptr) {
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1, nout2, nout3, nz4c});
  } else if (padm != nullptr) {
    AdvanceParticleRestartRealArrayOffset(data_size_, {nout1, nout2, nout3, nadm});
  }

  if (data_size_ != data_size) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "CC data size read from restart file not equal to size "
              << "of Hydro, MHD, Rad, and/or Z4c arrays, restart file is broken."
              << std::endl;
    restart_utils::AbortOnFatalError();
  }
  if (single_file_per_rank) {
    pm->restart_meta.common_prefix_bytes = static_cast<std::uint64_t>(headeroffset);
    pm->ValidateRestartShardCommonPrefix();
  }

  HostArray5D<Real> ccin("rst-cc-in", 1, 1, 1, 1, 1);
  HostFaceFld4D<Real> fcin("rst-fc-in", 1, 1, 1, 1);

  if (single_file_per_rank) {
    LoadSingleFileRestartData(pm, headeroffset, data_size_, nout1, nout2, nout3,
                              nhydro, nmhd, nrad, nforce, nz4c, nadm,
                              ccin, fcin);
  } else {
    // read CC data into host array
    int mygids = pm->gids_eachrank[global_variable::my_rank];
    if (mygids < 0) AbortParticleRestartLayoutOverflow();
    IOWrapperSizeT offset_myrank = headeroffset;
    if (!single_file_per_rank) {
      AdvanceParticleRestartOffset(offset_myrank, static_cast<IOWrapperSizeT>(mygids),
                                   data_size_);
    }
    IOWrapperSizeT myoffset = offset_myrank;

    // calculate max/min number of MeshBlocks across all ranks
    int noutmbs_max = pm->nmb_eachrank[0];
    int noutmbs_min = pm->nmb_eachrank[0];
    for (int i=0; i<(global_variable::nranks); ++i) {
      noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
      noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
    }

  if (phydro != nullptr) {
    ValidateParticleRestartAllocation({nmb, nhydro, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nhydro, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset, single_file_per_rank)
            != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC hydro data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset, single_file_per_rank)
            != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC hydro data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);
      }
    }
    Kokkos::deep_copy(Kokkos::subview(phydro->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nhydro});
    myoffset = offset_myrank;
  }

  if (pmhd != nullptr) {
    ValidateParticleRestartAllocation({nmb, nmhd, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nmhd, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                   Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset, single_file_per_rank)
            != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC mhd data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);
      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset, single_file_per_rank)
            != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC mhd data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pmhd->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nmhd});
    myoffset = offset_myrank;

    ValidateParticleRestartAllocation({nmb, nout3, nout2, nout1p1});
    Kokkos::realloc(fcin.x1f, nmb, nout3, nout2, nout1p1);
    ValidateParticleRestartAllocation({nmb, nout3, nout2p1, nout1});
    Kokkos::realloc(fcin.x2f, nmb, nout3, nout2p1, nout1);
    ValidateParticleRestartAllocation({nmb, nout3p1, nout2, nout1});
    Kokkos::realloc(fcin.x3f, nmb, nout3p1, nout2, nout1);
    // read FC data into host array, again one MeshBlock at a time
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(fcin.x1f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        auto fldcnt = x1fptr.size();

        if (resfile.Read_Reals_at_all(x1fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x1f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, static_cast<IOWrapperSizeT>(fldcnt),
                                     sizeof(Real));

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(fcin.x2f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = x2fptr.size();

        if (resfile.Read_Reals_at_all(x2fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x2f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, static_cast<IOWrapperSizeT>(fldcnt),
                                     sizeof(Real));

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(fcin.x3f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = x3fptr.size();

        if (resfile.Read_Reals_at_all(x3fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x3f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, static_cast<IOWrapperSizeT>(fldcnt),
                                     sizeof(Real));

        AdvanceParticleRestartPastFaceArrayRemainder(
            myoffset, data_size, {x1fptr.size(), x2fptr.size(), x3fptr.size()});
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(fcin.x1f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        auto fldcnt = x1fptr.size();

        if (resfile.Read_Reals_at(x1fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x1f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, static_cast<IOWrapperSizeT>(fldcnt),
                                     sizeof(Real));

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(fcin.x2f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = x2fptr.size();

        if (resfile.Read_Reals_at(x2fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x2f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, static_cast<IOWrapperSizeT>(fldcnt),
                                     sizeof(Real));

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(fcin.x3f, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        fldcnt = x3fptr.size();

        if (resfile.Read_Reals_at(x3fptr.data(), fldcnt, myoffset,
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Input b0.x3f field not read correctly from rst file, "
                << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, static_cast<IOWrapperSizeT>(fldcnt),
                                     sizeof(Real));

        AdvanceParticleRestartPastFaceArrayRemainder(
            myoffset, data_size, {x1fptr.size(), x2fptr.size(), x3fptr.size()});
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x1f);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x2f);
    Kokkos::deep_copy(Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL), fcin.x3f);
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1p1, nout2, nout3});
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1, nout2p1, nout3});
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3p1});
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    ValidateParticleRestartAllocation({nmb, nrad, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nrad, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC rad data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC rad data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);
      }
    }
    Kokkos::deep_copy(Kokkos::subview(prad->i0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nrad});
    myoffset = offset_myrank;
  }

  if (pturb != nullptr) {
    ValidateParticleRestartAllocation({nmb, nforce, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nforce, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC turb data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC turb data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pturb->force, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nforce});
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    ValidateParticleRestartAllocation({nmb, nz4c, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nz4c, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC z4c data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC z4c data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);
      }
    }
    Kokkos::deep_copy(Kokkos::subview(pz4c->u0, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nz4c});
    myoffset = offset_myrank;

    // We also need to reinitialize the ADM data.
    pz4c->Z4cToADM(pmy_mesh_->pmb_pack);
  } else if (padm != nullptr) {
    ValidateParticleRestartAllocation({nmb, nadm, nout3, nout2, nout1});
    Kokkos::realloc(ccin, nmb, nadm, nout3, nout2, nout1);
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to read, so read collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at_all(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC adm data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL);
        auto mbcnt = mbptr.size();
        if (resfile.Read_Reals_at(mbptr.data(), mbcnt, myoffset,
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "CC adm data not read correctly from rst file, "
                    << "restart file is broken." << std::endl;
          restart_utils::AbortOnFatalError();
        }
        AdvanceParticleRestartOffset(myoffset, 1, data_size);
      }
    }
    Kokkos::deep_copy(Kokkos::subview(padm->u_adm, std::make_pair(0,nmb), Kokkos::ALL,
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
    AdvanceParticleRestartRealArrayOffset(offset_myrank, {nout1, nout2, nout3, nadm});
    myoffset = offset_myrank;
  }
  }

  LoadParticleRestartData(pm, resfile, single_file_per_rank, headeroffset,
                          data_size_, nout1, nout2, nout3);

  // call problem generator again to re-initialize data, fn ptrs, as needed
#if USER_PROBLEM_ENABLED
  UserProblem(pin, true);
#else
  std::string pgen_fun_name = pin->GetOrAddString("problem", "pgen_name", "none");

  if (pgen_fun_name.compare("advection") == 0) {
    Advection(pin, true);
  } else if (pgen_fun_name.compare("cpaw") == 0) {
    AlfvenWave(pin, true);
  } else if (pgen_fun_name.compare("gr_bondi") == 0) {
    BondiAccretion(pin, true);
  } else if (pgen_fun_name.compare("tetrad") == 0) {
    CheckOrthonormalTetrad(pin, true);
  } else if (pgen_fun_name.compare("hohlraum") == 0) {
    Hohlraum(pin, true);
  } else if (pgen_fun_name.compare("linear_wave") == 0) {
    LinearWave(pin, true);
  } else if (pgen_fun_name.compare("implode") == 0) {
    LWImplode(pin, true);
  } else if (pgen_fun_name.compare("gr_monopole") == 0) {
    Monopole(pin, true);
  } else if (pgen_fun_name.compare("orszag_tang") == 0) {
    OrszagTang(pin, true);
  } else if (pgen_fun_name.compare("pic_paper_smooth_tsc_interface") == 0) {
    PICPaperSmoothTSCInterface(pin, true);
  } else if (pgen_fun_name.compare("pic_parallel_shock") == 0) {
    PICParallelShock(pin, true);
  } else if (pgen_fun_name.compare("q006_paper_multispecies_oscillation") == 0) {
    Q006PaperMultispeciesOscillation(pin, true);
  } else if (
      pgen_fun_name.compare("q006_paper_multispecies_oscillation_runtime_local")
      == 0) {
    Q006PaperMultispeciesOscillationRuntimeLocal(pin, true);
  } else if (pgen_fun_name.compare("q007_paper_crsi_linear_preparation") == 0) {
    Q007PaperCRSILinearPreparation(pin, true);
  } else if (pgen_fun_name.compare("q007_paper_crpai_linear_preparation") == 0) {
    Q007PaperCRPAILinearPreparation(pin, true);
  } else if (pgen_fun_name.compare("q023_paper_bell_linear") == 0) {
    Q023PaperBellLinear(pin, true);
  } else if (pgen_fun_name.compare("q043_bell_current_volume_aware") == 0) {
    Q043BellCurrentVolumeAware(pin, true);
  } else if (pgen_fun_name.compare("q029_hall_bell_linear") == 0) {
    Q029HallBellLinear(pin, true);
  } else if (pgen_fun_name.compare("q032_reduced_static_neutral_local") == 0) {
    Q032ReducedStaticNeutralLocal(pin, true);
  } else if (pgen_fun_name.compare("q033_crpai_transport_runtime_local") == 0) {
    Q033CRPAITransportRuntimeLocal(pin, true);
  } else if (pgen_fun_name.compare("rad_linear_wave") == 0) {
    RadiationLinearWave(pin, true);
  } else if (pgen_fun_name.compare("shock_tube") == 0) {
    ShockTube(pin, true);
  } else if (pgen_fun_name.compare("z4c_linear_wave") == 0) {
    Z4cLinearWave(pin, true);
  } else if (pgen_fun_name.compare("spherical_collapse") == 0) {
    SphericalCollapse(pin, true);
  } else if (pgen_fun_name.compare("diffusion") == 0) {
    Diffusion(pin, true);
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "Problem generator name could not be found in <problem> block in input file"
        << std::endl
        << "and it was not set by -D PROBLEM option on cmake command line during build"
        << std::endl
        << "Rerun cmake with -D PROBLEM=file to specify custom problem generator file"
        << std::endl;;
    restart_utils::AbortOnFatalError();
  }
#endif

  // Check that user defined BCs were enrolled if needed
  if (user_bcs) {
    if (user_bcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User BCs specified in <mesh> block, but not enrolled "
                << "during restart by SetProblemData()." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  // Check that user defined srcterms were enrolled if needed
  if (user_srcs) {
    if (user_srcs_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User SRCs specified in <problem> block, but not "
                << "enrolled by UserProblem()." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  // Check that user defined history outputs were enrolled if needed
  if (user_hist) {
    if (user_hist_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User history output specified in <problem> block, "
                << "but not enrolled by UserProblem()." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
  if (user_work_in_loop) {
    if (user_work_in_loop_func == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "User work-in-loop callback specified in <problem> "
                << "block, but not enrolled by UserProblem()." << std::endl;
      restart_utils::AbortOnFatalError();
    }
  }
}
