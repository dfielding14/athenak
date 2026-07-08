#ifndef BVALS_PARTICLE_COMPACTION_HPP_
#define BVALS_PARTICLE_COMPACTION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (see LICENSE file for details)
//========================================================================================
//! \file particle_compaction.hpp
//! \brief Host-side planning for particle migration and destruction compaction

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace particles {

// A receive replaces one of the first holes.  If fewer particles arrive than leave,
// surviving particles in [final_size, old_size) fill the remaining holes below
// final_size.  The two move ranges are disjoint, so the moves may execute in parallel.
struct ParticleCompactionPlan {
  int final_size = 0;
  std::vector<int> holes;
  std::vector<int> move_sources;
  std::vector<int> move_destinations;
};

inline bool BuildParticleCompactionPlan(
    const int old_size, const int receive_count, const std::vector<int> &send_indices,
    const std::vector<int> &destroy_indices, ParticleCompactionPlan *plan,
    std::string *error) {
  auto fail = [error](const std::string &message) {
    if (error != nullptr) *error = message;
    return false;
  };
  if (plan == nullptr) return fail("null compaction plan");
  *plan = ParticleCompactionPlan{};
  if (old_size < 0 || receive_count < 0) {
    return fail("negative particle or receive count");
  }

  plan->holes.reserve(send_indices.size() + destroy_indices.size());
  plan->holes.insert(plan->holes.end(), send_indices.begin(), send_indices.end());
  plan->holes.insert(plan->holes.end(), destroy_indices.begin(), destroy_indices.end());
  std::sort(plan->holes.begin(), plan->holes.end());
  for (std::size_t n = 0; n < plan->holes.size(); ++n) {
    if (plan->holes[n] < 0 || plan->holes[n] >= old_size) {
      return fail("particle removal index is outside the live array");
    }
    if (n > 0 && plan->holes[n] == plan->holes[n - 1]) {
      return fail("duplicate particle removal index");
    }
  }

  const std::int64_t final_size = static_cast<std::int64_t>(old_size) +
                                  static_cast<std::int64_t>(receive_count) -
                                  static_cast<std::int64_t>(plan->holes.size());
  if (final_size < 0 || final_size > std::numeric_limits<int>::max()) {
    return fail("particle compaction produces an invalid array size");
  }
  plan->final_size = static_cast<int>(final_size);

  // Receives consume the lowest holes.  Only unfilled holes that survive truncation
  // need a tail particle.  Every source is known not to be a removal index.
  const std::size_t first_unfilled = std::min(
      static_cast<std::size_t>(receive_count), plan->holes.size());
  for (std::size_t n = first_unfilled; n < plan->holes.size(); ++n) {
    if (plan->holes[n] >= plan->final_size) break;
    plan->move_destinations.push_back(plan->holes[n]);
  }

  auto next_hole = std::lower_bound(plan->holes.begin(), plan->holes.end(),
                                    plan->final_size);
  for (int source = plan->final_size; source < old_size; ++source) {
    if (next_hole != plan->holes.end() && *next_hole == source) {
      ++next_hole;
    } else {
      plan->move_sources.push_back(source);
    }
  }
  if (plan->move_sources.size() != plan->move_destinations.size()) {
    return fail("particle compaction source/destination count mismatch");
  }
  return true;
}

}  // namespace particles

#endif  // BVALS_PARTICLE_COMPACTION_HPP_
