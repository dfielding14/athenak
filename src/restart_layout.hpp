#ifndef RESTART_LAYOUT_HPP_
#define RESTART_LAYOUT_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart_layout.hpp
//! \brief Checked byte-layout helpers for restart serialization and reconstruction.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace restart_layout {

using Size = std::uint64_t;
using FailureHandler = void (*)(const std::string &);

inline Size CheckedAdd(Size left, Size right, FailureHandler fail,
                       const std::string &context) {
  if (right > std::numeric_limits<Size>::max() - left) {
    fail(context + " overflows.");
  }
  return left + right;
}

inline Size CheckedMultiply(Size left, Size right, FailureHandler fail,
                            const std::string &context) {
  if (left != 0 && right > std::numeric_limits<Size>::max()/left) {
    fail(context + " overflows.");
  }
  return left*right;
}

inline Size CheckedSubtract(Size left, Size right, FailureHandler fail,
                            const std::string &context) {
  if (right > left) {
    fail(context + " underflows.");
  }
  return left - right;
}

inline Size CheckedOffset(Size header_bytes, Size block_bytes, Size block_start,
                          FailureHandler fail, const std::string &context) {
  return CheckedAdd(header_bytes,
                    CheckedMultiply(block_bytes, block_start, fail, context),
                    fail, context);
}

inline Size CheckedPayloadBytes(Size header_bytes, Size block_bytes, Size blocks,
                                FailureHandler fail, const std::string &context) {
  return CheckedOffset(header_bytes, block_bytes, blocks, fail, context);
}

inline Size CheckedNonNegative(int value, FailureHandler fail,
                               const std::string &context) {
  if (value < 0) {
    fail(context + " must be non-negative.");
  }
  return static_cast<Size>(value);
}

inline Size CheckedPositive(int value, FailureHandler fail,
                            const std::string &context) {
  if (value <= 0) {
    fail(context + " must be positive.");
  }
  return static_cast<Size>(value);
}

inline Size CheckedExtentWithGhosts(int active_extent, int ghost_zones,
                                    bool expand, FailureHandler fail,
                                    const std::string &context) {
  Size extent = CheckedPositive(active_extent, fail, context + " active extent");
  if (!expand) {
    return 1;
  }
  Size ghosts = CheckedNonNegative(ghost_zones, fail, context + " ghost zones");
  return CheckedAdd(extent, CheckedMultiply(2, ghosts, fail, context), fail, context);
}

inline Size CheckedSizeT(std::size_t value, FailureHandler fail,
                         const std::string &context) {
  if constexpr (std::numeric_limits<std::size_t>::digits >
                std::numeric_limits<Size>::digits) {
    if (value > static_cast<std::size_t>(std::numeric_limits<Size>::max())) {
      fail(context + " exceeds restart count range.");
    }
  }
  return static_cast<Size>(value);
}

inline std::size_t CheckedMemorySize(Size value, FailureHandler fail,
                                     const std::string &context) {
  if constexpr (std::numeric_limits<std::size_t>::digits <
                std::numeric_limits<Size>::digits) {
    if (value > static_cast<Size>(std::numeric_limits<std::size_t>::max())) {
      fail(context + " does not fit in memory.");
    }
  }
  return static_cast<std::size_t>(value);
}

struct PayloadInput {
  Size nout1;
  Size nout2;
  Size nout3;
  Size nhydro;
  Size nmhd;
  Size nrad;
  Size nforce;
  Size nz4c;
  Size nadm;
  Size real_bytes;
};

struct PayloadLayout {
  Size hydro_bytes = 0;
  Size mhd_bytes = 0;
  Size mhd_x1f_bytes = 0;
  Size mhd_x2f_bytes = 0;
  Size mhd_x3f_bytes = 0;
  Size radiation_bytes = 0;
  Size forcing_bytes = 0;
  Size z4c_bytes = 0;
  Size adm_bytes = 0;
  Size data_bytes = 0;
  Size mhd_face_stride_remainder_bytes = 0;

  Size hydro_offset = 0;
  Size mhd_offset = 0;
  Size mhd_x1f_offset = 0;
  Size mhd_x2f_offset = 0;
  Size mhd_x3f_offset = 0;
  Size radiation_offset = 0;
  Size forcing_offset = 0;
  Size z4c_offset = 0;
  Size adm_offset = 0;

  static PayloadLayout Build(const PayloadInput &input, FailureHandler fail) {
    PayloadLayout layout;
    const auto field_bytes = [&](Size nx1, Size nx2, Size nx3, Size count,
                                 const std::string &context) {
      Size elements = CheckedMultiply(nx1, nx2, fail, context);
      elements = CheckedMultiply(elements, nx3, fail, context);
      elements = CheckedMultiply(elements, count, fail, context);
      return CheckedMultiply(elements, input.real_bytes, fail, context);
    };
    const auto append = [&](Size bytes, Size *offset, Size *total,
                            const std::string &context) {
      *offset = *total;
      *total = CheckedAdd(*total, bytes, fail, context);
    };

    layout.hydro_bytes = field_bytes(
        input.nout1, input.nout2, input.nout3, input.nhydro, "restart hydro bytes");
    append(layout.hydro_bytes, &layout.hydro_offset, &layout.data_bytes,
           "restart payload bytes");
    layout.mhd_bytes = field_bytes(
        input.nout1, input.nout2, input.nout3, input.nmhd, "restart MHD bytes");
    append(layout.mhd_bytes, &layout.mhd_offset, &layout.data_bytes,
           "restart payload bytes");
    if (input.nmhd > 0) {
      layout.mhd_x1f_bytes = field_bytes(
          CheckedAdd(input.nout1, 1, fail, "restart MHD x1 extent"),
          input.nout2, input.nout3, 1, "restart MHD x1-face bytes");
      layout.mhd_x2f_bytes = field_bytes(
          input.nout1, CheckedAdd(input.nout2, 1, fail, "restart MHD x2 extent"),
          input.nout3, 1, "restart MHD x2-face bytes");
      layout.mhd_x3f_bytes = field_bytes(
          input.nout1, input.nout2,
          CheckedAdd(input.nout3, 1, fail, "restart MHD x3 extent"),
          1, "restart MHD x3-face bytes");
    }
    append(layout.mhd_x1f_bytes, &layout.mhd_x1f_offset, &layout.data_bytes,
           "restart payload bytes");
    append(layout.mhd_x2f_bytes, &layout.mhd_x2f_offset, &layout.data_bytes,
           "restart payload bytes");
    append(layout.mhd_x3f_bytes, &layout.mhd_x3f_offset, &layout.data_bytes,
           "restart payload bytes");
    layout.radiation_bytes = field_bytes(
        input.nout1, input.nout2, input.nout3, input.nrad, "restart radiation bytes");
    append(layout.radiation_bytes, &layout.radiation_offset, &layout.data_bytes,
           "restart payload bytes");
    layout.forcing_bytes = field_bytes(
        input.nout1, input.nout2, input.nout3, input.nforce, "restart forcing bytes");
    append(layout.forcing_bytes, &layout.forcing_offset, &layout.data_bytes,
           "restart payload bytes");
    layout.z4c_bytes = field_bytes(
        input.nout1, input.nout2, input.nout3, input.nz4c, "restart Z4c bytes");
    append(layout.z4c_bytes, &layout.z4c_offset, &layout.data_bytes,
           "restart payload bytes");
    layout.adm_bytes = field_bytes(
        input.nout1, input.nout2, input.nout3, input.nadm, "restart ADM bytes");
    append(layout.adm_bytes, &layout.adm_offset, &layout.data_bytes,
           "restart payload bytes");
    Size mhd_face_bytes = CheckedAdd(layout.mhd_x1f_bytes, layout.mhd_x2f_bytes,
                                     fail, "restart MHD face bytes");
    mhd_face_bytes = CheckedAdd(mhd_face_bytes, layout.mhd_x3f_bytes, fail,
                                "restart MHD face bytes");
    layout.mhd_face_stride_remainder_bytes = CheckedSubtract(
        layout.data_bytes, mhd_face_bytes, fail, "restart MHD face stride remainder");
    return layout;
  }
};

struct ManifestBudget {
  Size bytes = 0;
  Size max_bytes;

  void Add(Size record_bytes, FailureHandler fail) {
    bytes = CheckedAdd(bytes, record_bytes, fail, "node restart manifest bytes");
    if (bytes > max_bytes) {
      fail("node restart manifest exceeds the " + std::to_string(max_bytes) +
           "-byte limit.");
    }
  }
};

struct HeaderInput {
  Size parameter_bytes;
  Size marker_bytes;
  Size fixed_mesh_bytes;
  Size meshblocks;
  Size metadata_record_bytes;
  Size state_bytes;
  Size data_size_record_bytes;
};

struct HeaderLayout {
  Size metadata_bytes = 0;
  Size header_bytes = 0;

  static HeaderLayout Build(const HeaderInput &input, FailureHandler fail) {
    HeaderLayout layout;
    layout.metadata_bytes = CheckedMultiply(
        input.meshblocks, input.metadata_record_bytes, fail,
        "restart MeshBlock metadata bytes");
    layout.header_bytes = CheckedAdd(input.parameter_bytes, input.marker_bytes, fail,
                                     "restart header bytes");
    layout.header_bytes = CheckedAdd(layout.header_bytes, input.fixed_mesh_bytes, fail,
                                     "restart header bytes");
    layout.header_bytes = CheckedAdd(layout.header_bytes, layout.metadata_bytes, fail,
                                     "restart header bytes");
    layout.header_bytes = CheckedAdd(layout.header_bytes, input.state_bytes, fail,
                                     "restart header bytes");
    layout.header_bytes = CheckedAdd(layout.header_bytes, input.data_size_record_bytes,
                                     fail, "restart header bytes");
    return layout;
  }
};

}  // namespace restart_layout

#endif  // RESTART_LAYOUT_HPP_
