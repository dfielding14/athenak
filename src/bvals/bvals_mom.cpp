//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file bvals_mom.cpp
//! \brief Host transport for duplicated paper_smooth AMR particle moment records

#include <algorithm>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "bvals.hpp"

namespace particles {
namespace {
using Record = PaperSmoothMomentRecord;

static_assert(std::is_standard_layout<Record>::value,
              "paper_smooth moment record must have standard layout");
static_assert(std::is_trivial<Record>::value,
              "paper_smooth moment record must remain trivial");
static_assert(std::is_trivially_copyable<Record>::value,
              "paper_smooth moment record must be byte-copyable");
static_assert(sizeof(Record) == 4*sizeof(std::uint32_t) + 12*sizeof(Real),
              "paper_smooth moment record must not contain implicit padding");

#if MPI_PARALLEL_ENABLED
bool BuildDisplacements(const std::vector<int> &counts, std::vector<int> &displs,
                        int &total) {
  total = 0;
  displs.resize(counts.size());
  for (std::size_t n = 0; n < counts.size(); ++n) {
    if (counts[n] < 0 || counts[n] > (std::numeric_limits<int>::max() - total)) {
      return false;
    }
    displs[n] = total;
    total += counts[n];
  }
  return true;
}

bool AllRanksValid(const bool local_valid, MPI_Comm comm) {
  int local = local_valid ? 1 : 0;
  int global = 0;
  return (MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, comm) == MPI_SUCCESS &&
          global != 0);
}
#endif
} // namespace

//----------------------------------------------------------------------------------------
//! \brief Create an isolated communicator for the synchronous paper_smooth exchange.

PaperSmoothMomentRecordTransport::PaperSmoothMomentRecordTransport() :
    queued_key_count_(0),
    delivered_key_count_(0),
#if MPI_PARALLEL_ENABLED
    mpi_comm_mom_(MPI_COMM_NULL),
    mpi_ready_(false),
#endif
    my_rank_(0),
    nranks_(1) {
#if MPI_PARALLEL_ENABLED
  int initialized = 0;
  int finalized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized == 0) return;
  if (MPI_Finalized(&finalized) != MPI_SUCCESS || finalized != 0) return;
  if (MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm_mom_) != MPI_SUCCESS) return;
  if (MPI_Comm_set_errhandler(mpi_comm_mom_, MPI_ERRORS_RETURN) != MPI_SUCCESS) return;
  if (MPI_Comm_rank(mpi_comm_mom_, &my_rank_) != MPI_SUCCESS) return;
  if (MPI_Comm_size(mpi_comm_mom_, &nranks_) != MPI_SUCCESS) return;
  mpi_ready_ = true;
#endif
}

//----------------------------------------------------------------------------------------
//! \brief Free the isolated communicator when MPI is still active.

PaperSmoothMomentRecordTransport::~PaperSmoothMomentRecordTransport() {
#if MPI_PARALLEL_ENABLED
  int initialized = 0;
  int finalized = 0;
  if (mpi_comm_mom_ == MPI_COMM_NULL) return;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized == 0) return;
  if (MPI_Finalized(&finalized) != MPI_SUCCESS || finalized != 0) return;
  MPI_Comm_free(&mpi_comm_mom_);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief Clear records queued or delivered during the previous deposition stage.

void PaperSmoothMomentRecordTransport::Reset() {
  records_.clear();
  outgoing_.clear();
  send_records_.clear();
  recv_records_.clear();
  ResetKeys(queued_keys_, queued_key_count_);
  ResetKeys(delivered_keys_, delivered_key_count_);
}

//----------------------------------------------------------------------------------------
//! \brief Queue one unique cross-level duplicate without touching migration data.

PaperSmoothRecordStatus PaperSmoothMomentRecordTransport::Queue(
    const PaperSmoothMomentRecord &record, const int dest_rank,
    const int source_level, const int dest_level) {
  if (source_level == dest_level) return PaperSmoothRecordStatus::duplicate;
  return QueueReceiver(record, dest_rank);
}

//----------------------------------------------------------------------------------------
//! \brief Queue one unique receiver-leaf record without touching migration data.

PaperSmoothRecordStatus PaperSmoothMomentRecordTransport::QueueReceiver(
    const PaperSmoothMomentRecord &record, const int dest_rank) {
  if (record.dest_gid < 0 || record.ptag < 0 ||
      !PaperSmoothImageCodeValid(record.reserved)) {
    return PaperSmoothRecordStatus::invalid;
  }
  if (dest_rank < 0 || dest_rank >= nranks_) return PaperSmoothRecordStatus::invalid;
  if (!InsertKey(queued_keys_, queued_key_count_, MakeKey(record))) {
    return PaperSmoothRecordStatus::duplicate;
  }

  if (dest_rank == my_rank_) return AddRecord(record);
  outgoing_.push_back({dest_rank, record});
  return PaperSmoothRecordStatus::accepted;
}

//----------------------------------------------------------------------------------------
//! \brief Retain local records and synchronously exchange records bound for remote ranks.

bool PaperSmoothMomentRecordTransport::Exchange() {
#if MPI_PARALLEL_ENABLED
  if (!mpi_ready_) return false;

  constexpr int record_size = static_cast<int>(sizeof(PaperSmoothMomentRecord));
  send_counts_.assign(nranks_, 0);
  recv_counts_.assign(nranks_, 0);
  bool local_valid = true;
  for (const auto &entry : outgoing_) {
    if (entry.dest_rank < 0 || entry.dest_rank >= nranks_ ||
        send_counts_[entry.dest_rank] > (std::numeric_limits<int>::max() - record_size)) {
      local_valid = false;
      break;
    }
    send_counts_[entry.dest_rank] += record_size;
  }

  int send_bytes = 0;
  local_valid = local_valid && BuildDisplacements(send_counts_, send_displs_, send_bytes);
  if (!AllRanksValid(local_valid, mpi_comm_mom_)) return false;

  send_records_.resize(static_cast<std::size_t>(send_bytes/record_size));
  std::vector<std::size_t> offsets(static_cast<std::size_t>(nranks_));
  for (int rank = 0; rank < nranks_; ++rank) {
    offsets[rank] = static_cast<std::size_t>(send_displs_[rank]/record_size);
  }
  for (const auto &entry : outgoing_) {
    send_records_[offsets[entry.dest_rank]++] = entry.record;
  }

  if (MPI_Alltoall(send_counts_.data(), 1, MPI_INT, recv_counts_.data(), 1, MPI_INT,
                   mpi_comm_mom_) != MPI_SUCCESS) {
    return false;
  }

  int recv_bytes = 0;
  local_valid = BuildDisplacements(recv_counts_, recv_displs_, recv_bytes);
  for (const int count : recv_counts_) {
    if ((count % record_size) != 0) local_valid = false;
  }
  if (!AllRanksValid(local_valid, mpi_comm_mom_)) return false;

  recv_records_.resize(static_cast<std::size_t>(recv_bytes/record_size));
  if (MPI_Alltoallv(send_records_.data(), send_counts_.data(), send_displs_.data(),
                    MPI_BYTE, recv_records_.data(), recv_counts_.data(),
                    recv_displs_.data(), MPI_BYTE, mpi_comm_mom_) != MPI_SUCCESS) {
    return false;
  }

  outgoing_.clear();
  local_valid = true;
  for (const auto &record : recv_records_) {
    if (AddRecord(record) == PaperSmoothRecordStatus::invalid) local_valid = false;
  }
  if (!AllRanksValid(local_valid, mpi_comm_mom_)) return false;
#else
  outgoing_.clear();
#endif
  return true;
}

//----------------------------------------------------------------------------------------
//! \brief Return retained bytes used by vectors that store records.

std::uint64_t PaperSmoothMomentRecordTransport::RecordAllocationBytes() const {
  return static_cast<std::uint64_t>(records_.capacity()) *
             sizeof(PaperSmoothMomentRecord) +
         static_cast<std::uint64_t>(outgoing_.capacity()) * sizeof(DestinationRecord) +
         static_cast<std::uint64_t>(send_records_.capacity() + recv_records_.capacity()) *
             sizeof(PaperSmoothMomentRecord);
}

//----------------------------------------------------------------------------------------
//! \brief Return retained bytes used by deduplication and MPI metadata vectors.

std::uint64_t PaperSmoothMomentRecordTransport::MetadataAllocationBytes() const {
  std::uint64_t bytes =
      static_cast<std::uint64_t>(queued_keys_.capacity() + delivered_keys_.capacity()) *
      sizeof(RecordKey);
#if MPI_PARALLEL_ENABLED
  bytes += static_cast<std::uint64_t>(send_counts_.capacity() + recv_counts_.capacity() +
                                     send_displs_.capacity() + recv_displs_.capacity()) *
           sizeof(int);
#endif
  return bytes;
}

//----------------------------------------------------------------------------------------
//! \brief Return retained host payload allocation bytes owned by this helper.

std::uint64_t PaperSmoothMomentRecordTransport::AllocationBytes() const {
  return RecordAllocationBytes() + MetadataAllocationBytes();
}

//----------------------------------------------------------------------------------------
//! \brief Pack the stable particle provenance key and destination GID for deduplication.

PaperSmoothMomentRecordTransport::RecordKey
PaperSmoothMomentRecordTransport::MakeKey(const PaperSmoothMomentRecord &record) {
  return {static_cast<std::uint32_t>(record.ptag),
          static_cast<std::uint32_t>(record.dest_gid), record.reserved};
}

//----------------------------------------------------------------------------------------
//! \brief Hash a packed record key using the SplitMix64 finalizer.

std::uint64_t PaperSmoothMomentRecordTransport::HashKey(const RecordKey &record_key) {
  std::uint64_t key = record_key.ptag;
  key ^= static_cast<std::uint64_t>(record_key.dest_gid) +
         UINT64_C(0x9e3779b97f4a7c15) + (key << 6) + (key >> 2);
  key ^= static_cast<std::uint64_t>(record_key.image_code) +
         UINT64_C(0x9e3779b97f4a7c15) + (key << 6) + (key >> 2);
  key = (key ^ (key >> 30))*UINT64_C(0xbf58476d1ce4e5b9);
  key = (key ^ (key >> 27))*UINT64_C(0x94d049bb133111eb);
  return key ^ (key >> 31);
}

//----------------------------------------------------------------------------------------
//! \brief Compare stable particle, receiver, and periodic-image record identities.

bool PaperSmoothMomentRecordTransport::KeysEqual(const RecordKey &lhs,
                                                 const RecordKey &rhs) {
  return lhs.ptag == rhs.ptag && lhs.dest_gid == rhs.dest_gid &&
         lhs.image_code == rhs.image_code;
}

//----------------------------------------------------------------------------------------
//! \brief Return true for an unused open-addressing table slot.

bool PaperSmoothMomentRecordTransport::IsEmptyKey(const RecordKey &key) {
  return key.ptag == empty_key_component_;
}

//----------------------------------------------------------------------------------------
//! \brief Resize an open-addressed key table while preserving existing keys.

void PaperSmoothMomentRecordTransport::RehashKeys(std::vector<RecordKey> &table,
                                                  const std::size_t new_size) {
  const RecordKey empty_key{empty_key_component_, empty_key_component_,
                            empty_key_component_};
  std::vector<RecordKey> new_table(new_size, empty_key);
  for (const RecordKey &key : table) {
    if (IsEmptyKey(key)) continue;
    std::size_t slot = static_cast<std::size_t>(HashKey(key)) & (new_size - 1);
    while (!IsEmptyKey(new_table[slot])) slot = (slot + 1) & (new_size - 1);
    new_table[slot] = key;
  }
  table.swap(new_table);
}

//----------------------------------------------------------------------------------------
//! \brief Insert a key into an allocation-accountable open-addressed deduplication table.

bool PaperSmoothMomentRecordTransport::InsertKey(std::vector<RecordKey> &table,
                                                 std::size_t &count,
                                                 const RecordKey &key) {
  if (table.empty()) {
    RehashKeys(table, 16);
  } else if ((count + 1)*10 > table.size()*7) {
    RehashKeys(table, 2*table.size());
  }

  std::size_t slot = static_cast<std::size_t>(HashKey(key)) & (table.size() - 1);
  while (!IsEmptyKey(table[slot])) {
    if (KeysEqual(table[slot], key)) return false;
    slot = (slot + 1) & (table.size() - 1);
  }
  table[slot] = key;
  ++count;
  return true;
}

//----------------------------------------------------------------------------------------
//! \brief Reset a key table while retaining its allocation for the next stage.

void PaperSmoothMomentRecordTransport::ResetKeys(std::vector<RecordKey> &table,
                                                 std::size_t &count) {
  const RecordKey empty_key{empty_key_component_, empty_key_component_,
                            empty_key_component_};
  std::fill(table.begin(), table.end(), empty_key);
  count = 0;
}

//----------------------------------------------------------------------------------------
//! \brief Retain a local or received record once per provenance key and destination GID.

PaperSmoothRecordStatus PaperSmoothMomentRecordTransport::AddRecord(
    const PaperSmoothMomentRecord &record) {
  if (record.dest_gid < 0 || record.ptag < 0 ||
      !PaperSmoothImageCodeValid(record.reserved)) {
    return PaperSmoothRecordStatus::invalid;
  }
  if (!InsertKey(delivered_keys_, delivered_key_count_, MakeKey(record))) {
    return PaperSmoothRecordStatus::duplicate;
  }
  records_.push_back(record);
  return PaperSmoothRecordStatus::accepted;
}
} // namespace particles
