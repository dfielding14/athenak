#ifndef PARTICLES_PARTICLES_DATA_STRUCTS_HPP_
#define PARTICLES_PARTICLES_DATA_STRUCTS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#include <cstdint>

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \struct ParticleLocationData
//! \brief data describing location of data for particles communicated with MPI

struct ParticleLocationData {
  int prtcl_indx;   // index in particle array
  int dest_gid;     // GID of target MeshBlock
  int dest_rank;    // rank of target MeshBlock
};

// Custom operators to sort ParticleLocationData array by dest_rank or prtcl_indx
struct {
  bool operator()(ParticleLocationData a, ParticleLocationData b)
    const { return a.dest_rank < b.dest_rank; }
} SortByRank;
struct {
  bool operator()(ParticleLocationData a, ParticleLocationData b)
    const { return a.prtcl_indx < b.prtcl_indx; }
} SortByIndex;

//----------------------------------------------------------------------------------------
//! \struct ParticleMessageData
//! \brief Data describing MPI messages containing particles

struct ParticleMessageData {
  int sendrank;  // rank of sender
  int recvrank;  // rank of receiver
  int nprtcls;   // number of particles in message
  ParticleMessageData(int a, int b, int c) :
    sendrank(a), recvrank(b), nprtcls(c) {}
};

namespace particles {
//----------------------------------------------------------------------------------------
//! \struct PaperSmoothMomentRecord
//! \brief Host-transport record for receiver-resolution paper_smooth AMR deposition

constexpr std::uint32_t kPaperSmoothDepositRhoJ = 1U << 0;
constexpr std::uint32_t kPaperSmoothDepositEBDot = 1U << 1;
constexpr std::uint32_t kPaperSmoothDepositMomentumFeedback = 1U << 2;
constexpr std::uint32_t kPaperSmoothDepositEnergyFeedback = 1U << 3;
constexpr std::uint32_t kPaperSmoothImageAxisMask = 3U;
constexpr std::uint32_t kPaperSmoothImageCodeMask = 0x3fU;

KOKKOS_INLINE_FUNCTION
constexpr std::uint32_t PaperSmoothEncodeImageShift(const int shift) {
  return (shift == -1) ? 2U : ((shift == 1) ? 1U : ((shift == 0) ? 0U : 3U));
}

KOKKOS_INLINE_FUNCTION
constexpr std::uint32_t PaperSmoothEncodeImageCode(const int shift_x1,
                                                   const int shift_x2,
                                                   const int shift_x3) {
  return PaperSmoothEncodeImageShift(shift_x1) |
         (PaperSmoothEncodeImageShift(shift_x2) << 2) |
         (PaperSmoothEncodeImageShift(shift_x3) << 4);
}

KOKKOS_INLINE_FUNCTION
constexpr int PaperSmoothDecodeImageShift(const std::uint32_t code,
                                          const int axis) {
  const std::uint32_t encoded = (code >> (2*axis)) & kPaperSmoothImageAxisMask;
  return (encoded == 2U) ? -1 : ((encoded == 1U) ? 1 : 0);
}

KOKKOS_INLINE_FUNCTION
constexpr bool PaperSmoothImageCodeValid(const std::uint32_t code) {
  return ((code & ~kPaperSmoothImageCodeMask) == 0U) &&
         ((code & kPaperSmoothImageAxisMask) != 3U) &&
         (((code >> 2) & kPaperSmoothImageAxisMask) != 3U) &&
         (((code >> 4) & kPaperSmoothImageAxisMask) != 3U);
}

struct PaperSmoothMomentRecord {
  std::int32_t dest_gid;
  std::int32_t ptag;
  std::int32_t source_rank;
  std::int32_t source_index; // unique with source_rank during one deposition stage
  std::uint32_t deposit_flags;
  std::uint32_t reserved; // periodic receiver-image code; zero is the native image
  Real x, y, z;
  Real q_macro; // extensive macro-charge
  Real vx, vy, vz;
  Real ebdot;
  Real dpxdt, dpydt, dpzdt, dedt; // extensive feedback rates
};
} // namespace particles

#endif // PARTICLES_PARTICLES_DATA_STRUCTS_HPP_
