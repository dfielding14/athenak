//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file io_wrapper.cpp
//! \brief functions that provide wrapper for MPI-IO versus serial input/output

#include "io_wrapper.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"

namespace {

constexpr const char* kTestMaxMpiBytesEnv = "ATHENAK_TEST_MAX_MPI_BYTES";

[[noreturn]] void FatalIOError(const std::string& message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
            << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

IOWrapperSizeT CheckedByteCount(IOWrapperSizeT size, IOWrapperSizeT count,
                                const char* context) {
  if (size != 0 && count > std::numeric_limits<IOWrapperSizeT>::max()/size) {
    FatalIOError(std::string(context) + " byte count overflow.");
  }
  return size*count;
}

IOWrapperSizeT CheckedOffsetAdd(IOWrapperSizeT offset, IOWrapperSizeT increment,
                                const char* context) {
  if (increment > std::numeric_limits<IOWrapperSizeT>::max() - offset) {
    FatalIOError(std::string(context) + " offset overflow.");
  }
  return offset + increment;
}

std::size_t CheckedSizeT(IOWrapperSizeT value, const char* context) {
  if (value > std::numeric_limits<std::size_t>::max()) {
    FatalIOError(std::string(context) + " exceeds std::size_t range.");
  }
  return static_cast<std::size_t>(value);
}

std::size_t CompletedElements(IOWrapperSizeT bytes, IOWrapperSizeT element_size,
                              const char* context) {
  return CheckedSizeT(bytes/element_size, context);
}

std::size_t SizeOfDatatype(const std::string& datatype) {
  if (datatype == "byte") {
    return sizeof(char);
  }
  if (datatype == "int") {
    return sizeof(int);
  }
  if (datatype == "float") {
    return sizeof(float);
  }
  if (datatype == "double") {
    return sizeof(double);
  }
  if (datatype == "Real") {
    return sizeof(Real);
  }
  FatalIOError("Unrecognized datatype '" + datatype + "'.");
}

#if MPI_PARALLEL_ENABLED
std::string MpiErrorMessage(const char* context, int mpi_error) {
  char message[MPI_MAX_ERROR_STRING];
  int message_length = 0;
  MPI_Error_string(mpi_error, message, &message_length);
  return std::string(context) + " failed with MPI error: " +
         std::string(message, message_length);
}

[[noreturn]] void FatalMpiError(const char* context, int mpi_error) {
  FatalIOError(MpiErrorMessage(context, mpi_error));
}

void PrintMpiError(const char* context, int mpi_error) {
  std::cerr << MpiErrorMessage(context, mpi_error) << std::endl;
}

IOWrapperSizeT GetConfiguredMaxMpiBytes() {
  static const IOWrapperSizeT max_chunk_bytes = []() {
    const IOWrapperSizeT mpi_max =
        static_cast<IOWrapperSizeT>(std::numeric_limits<int>::max());
    const char* value = std::getenv(kTestMaxMpiBytesEnv);
    if (value == nullptr) {
      return mpi_max;
    }
    if (*value == '\0') {
      FatalIOError(std::string(kTestMaxMpiBytesEnv) +
                   " must be a decimal integer between 1 and INT_MAX.");
    }
    for (const char* digit = value; *digit != '\0'; ++digit) {
      if (*digit < '0' || *digit > '9') {
        FatalIOError(std::string(kTestMaxMpiBytesEnv) +
                     " must be a decimal integer between 1 and INT_MAX.");
      }
    }
    errno = 0;
    char* parse_end = nullptr;
    std::uint64_t parsed = std::strtoull(value, &parse_end, 10);
    if (errno == ERANGE || parse_end == value || *parse_end != '\0' ||
        parsed == 0 || parsed > mpi_max) {
      FatalIOError(std::string(kTestMaxMpiBytesEnv) +
                   " must be a decimal integer between 1 and INT_MAX.");
    }
    return static_cast<IOWrapperSizeT>(parsed);
  }();
  return max_chunk_bytes;
}

MPI_Offset CheckedMpiOffset(IOWrapperSizeT offset, const char* context) {
  if constexpr (std::numeric_limits<MPI_Offset>::digits <
                std::numeric_limits<IOWrapperSizeT>::digits) {
    if (offset > static_cast<IOWrapperSizeT>(
                     std::numeric_limits<MPI_Offset>::max())) {
      FatalIOError(std::string(context) + " exceeds MPI_Offset range.");
    }
  }
  return static_cast<MPI_Offset>(offset);
}

void PreflightPositionedMpiRange(IOWrapperSizeT offset,
                                 IOWrapperSizeT total_bytes,
                                 const char* context) {
  if (total_bytes == 0) {
    return;
  }
  IOWrapperSizeT end_offset =
      CheckedOffsetAdd(offset, total_bytes - 1, context);
  CheckedMpiOffset(end_offset, context);
}

IOWrapperSizeT TransferredBytes(MPI_Status* status, const char* context,
                                bool* valid) {
  int transferred = 0;
  int mpi_error = MPI_Get_count(status, MPI_BYTE, &transferred);
  if (mpi_error != MPI_SUCCESS) {
    PrintMpiError(context, mpi_error);
    *valid = false;
    return 0;
  }
  if (transferred == MPI_UNDEFINED || transferred < 0) {
    std::cerr << context << " returned an invalid byte count." << std::endl;
    *valid = false;
    return 0;
  }
  *valid = true;
  return static_cast<IOWrapperSizeT>(transferred);
}

bool CollectiveRoundSucceeded(MPI_Comm comm, bool local_success,
                              const char* context) {
  int local_failed = local_success ? 0 : 1;
  int any_failed = 0;
  int mpi_error = MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX,
                                comm);
  if (mpi_error != MPI_SUCCESS) {
    FatalMpiError(context, mpi_error);
  }
  return any_failed == 0;
}

IOWrapperSizeT CollectiveMaxBytes(MPI_Comm comm, IOWrapperSizeT local_bytes,
                                  const char* context) {
  IOWrapperSizeT max_bytes = 0;
  int mpi_error = MPI_Allreduce(&local_bytes, &max_bytes, 1, MPI_UINT64_T,
                                MPI_MAX, comm);
  if (mpi_error != MPI_SUCCESS) {
    FatalMpiError(context, mpi_error);
  }
  return max_bytes;
}

IOWrapperSizeT AgreedMaxMpiBytes(MPI_Comm comm, const char* context) {
  IOWrapperSizeT local_bytes = GetConfiguredMaxMpiBytes();
  IOWrapperSizeT min_bytes = 0;
  IOWrapperSizeT max_bytes = 0;
  int mpi_error = MPI_Allreduce(&local_bytes, &min_bytes, 1, MPI_UINT64_T,
                                MPI_MIN, comm);
  if (mpi_error != MPI_SUCCESS) {
    FatalMpiError(context, mpi_error);
  }
  mpi_error = MPI_Allreduce(&local_bytes, &max_bytes, 1, MPI_UINT64_T, MPI_MAX,
                            comm);
  if (mpi_error != MPI_SUCCESS) {
    FatalMpiError(context, mpi_error);
  }
  if (min_bytes != max_bytes) {
    FatalIOError(std::string(kTestMaxMpiBytesEnv) +
                 " must have the same value on every communicator rank.");
  }
  return min_bytes;
}

IOWrapperSizeT ChunkedMpiByteRead(IOWrapperFile fh, char* buf,
                                  IOWrapperSizeT total_bytes) {
  IOWrapperSizeT bytes_read = 0;
  IOWrapperSizeT max_chunk = GetConfiguredMaxMpiBytes();
  while (bytes_read < total_bytes) {
    IOWrapperSizeT chunk_bytes = std::min(max_chunk, total_bytes - bytes_read);
    MPI_Status status;
    int mpi_error = MPI_File_read(
        fh, buf + CheckedSizeT(bytes_read, "Read_bytes buffer offset"),
        static_cast<int>(chunk_bytes), MPI_BYTE, &status);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError("MPI_File_read", mpi_error);
      return 0;
    }
    bool valid = false;
    IOWrapperSizeT transferred = TransferredBytes(
        &status, "MPI_File_read MPI_Get_count", &valid);
    if (!valid) {
      return 0;
    }
    bytes_read = CheckedOffsetAdd(bytes_read, transferred, "MPI_File_read");
    if (transferred != chunk_bytes) {
      break;
    }
  }
  return bytes_read;
}

IOWrapperSizeT ChunkedMpiByteReadAt(IOWrapperFile fh, char* buf,
                                    IOWrapperSizeT total_bytes,
                                    IOWrapperSizeT offset) {
  PreflightPositionedMpiRange(offset, total_bytes, "MPI_File_read_at");
  IOWrapperSizeT bytes_read = 0;
  IOWrapperSizeT max_chunk = GetConfiguredMaxMpiBytes();
  while (bytes_read < total_bytes) {
    IOWrapperSizeT chunk_bytes = std::min(max_chunk, total_bytes - bytes_read);
    MPI_Offset mpi_offset = CheckedMpiOffset(
        CheckedOffsetAdd(offset, bytes_read, "MPI_File_read_at"),
        "MPI_File_read_at");
    MPI_Status status;
    int mpi_error = MPI_File_read_at(
        fh, mpi_offset,
        buf + CheckedSizeT(bytes_read, "Read_bytes_at buffer offset"),
        static_cast<int>(chunk_bytes), MPI_BYTE, &status);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError("MPI_File_read_at", mpi_error);
      return 0;
    }
    bool valid = false;
    IOWrapperSizeT transferred = TransferredBytes(
        &status, "MPI_File_read_at MPI_Get_count", &valid);
    if (!valid) {
      return 0;
    }
    bytes_read = CheckedOffsetAdd(bytes_read, transferred, "MPI_File_read_at");
    if (transferred != chunk_bytes) {
      break;
    }
  }
  return bytes_read;
}

IOWrapperSizeT ChunkedMpiByteReadAtAll(IOWrapperFile fh, MPI_Comm comm,
                                       char* buf, IOWrapperSizeT total_bytes,
                                       IOWrapperSizeT offset) {
  PreflightPositionedMpiRange(offset, total_bytes, "MPI_File_read_at_all");
  IOWrapperSizeT max_total_bytes =
      CollectiveMaxBytes(comm, total_bytes, "MPI_Allreduce read byte count");
  IOWrapperSizeT max_chunk =
      AgreedMaxMpiBytes(comm, "MPI_Allreduce read chunk limit");
  char dummy = '\0';
  IOWrapperSizeT bytes_read = 0;
  IOWrapperSizeT chunk_begin = 0;
  while (chunk_begin < max_total_bytes) {
    IOWrapperSizeT scheduled_bytes =
        std::min(max_chunk, max_total_bytes - chunk_begin);
    IOWrapperSizeT local_bytes = 0;
    char* local_buf = &dummy;
    MPI_Offset mpi_offset = 0;
    if (chunk_begin < total_bytes) {
      local_bytes = std::min(scheduled_bytes, total_bytes - chunk_begin);
      local_buf = buf + CheckedSizeT(chunk_begin,
                                     "Read_bytes_at_all buffer offset");
      mpi_offset = CheckedMpiOffset(
          CheckedOffsetAdd(offset, chunk_begin, "MPI_File_read_at_all"),
          "MPI_File_read_at_all");
    }
    MPI_Status status;
    int mpi_error = MPI_File_read_at_all(
        fh, mpi_offset, local_buf, static_cast<int>(local_bytes), MPI_BYTE,
        &status);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError("MPI_File_read_at_all", mpi_error);
    }
    bool valid = false;
    IOWrapperSizeT transferred = 0;
    if (mpi_error == MPI_SUCCESS) {
      transferred = TransferredBytes(
          &status, "MPI_File_read_at_all MPI_Get_count", &valid);
    }
    bool local_success =
        mpi_error == MPI_SUCCESS && valid && transferred == local_bytes;
    if (!CollectiveRoundSucceeded(comm, local_success,
                                  "MPI_Allreduce read result")) {
      return 0;
    }
    bytes_read = CheckedOffsetAdd(bytes_read, transferred,
                                  "MPI_File_read_at_all");
    chunk_begin = CheckedOffsetAdd(chunk_begin, scheduled_bytes,
                                   "MPI_File_read_at_all progress");
  }
  return bytes_read;
}

IOWrapperSizeT ChunkedMpiByteWrite(IOWrapperFile fh, const char* buf,
                                   IOWrapperSizeT total_bytes) {
  IOWrapperSizeT bytes_written = 0;
  IOWrapperSizeT max_chunk = GetConfiguredMaxMpiBytes();
  while (bytes_written < total_bytes) {
    IOWrapperSizeT chunk_bytes =
        std::min(max_chunk, total_bytes - bytes_written);
    MPI_Status status;
    int mpi_error = MPI_File_write(
        fh, const_cast<char*>(
                buf + CheckedSizeT(bytes_written, "Write_any_type buffer offset")),
        static_cast<int>(chunk_bytes), MPI_BYTE, &status);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError("MPI_File_write", mpi_error);
      return 0;
    }
    bool valid = false;
    IOWrapperSizeT transferred = TransferredBytes(
        &status, "MPI_File_write MPI_Get_count", &valid);
    if (!valid) {
      return 0;
    }
    bytes_written = CheckedOffsetAdd(bytes_written, transferred,
                                     "MPI_File_write");
    if (transferred != chunk_bytes) {
      break;
    }
  }
  return bytes_written;
}

IOWrapperSizeT ChunkedMpiByteWriteAt(IOWrapperFile fh, const char* buf,
                                     IOWrapperSizeT total_bytes,
                                     IOWrapperSizeT offset) {
  PreflightPositionedMpiRange(offset, total_bytes, "MPI_File_write_at");
  IOWrapperSizeT bytes_written = 0;
  IOWrapperSizeT max_chunk = GetConfiguredMaxMpiBytes();
  while (bytes_written < total_bytes) {
    IOWrapperSizeT chunk_bytes =
        std::min(max_chunk, total_bytes - bytes_written);
    MPI_Offset mpi_offset = CheckedMpiOffset(
        CheckedOffsetAdd(offset, bytes_written, "MPI_File_write_at"),
        "MPI_File_write_at");
    MPI_Status status;
    int mpi_error = MPI_File_write_at(
        fh, mpi_offset,
        const_cast<char*>(
            buf + CheckedSizeT(bytes_written, "Write_any_type_at buffer offset")),
        static_cast<int>(chunk_bytes), MPI_BYTE, &status);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError("MPI_File_write_at", mpi_error);
      return 0;
    }
    bool valid = false;
    IOWrapperSizeT transferred = TransferredBytes(
        &status, "MPI_File_write_at MPI_Get_count", &valid);
    if (!valid) {
      return 0;
    }
    bytes_written = CheckedOffsetAdd(bytes_written, transferred,
                                     "MPI_File_write_at");
    if (transferred != chunk_bytes) {
      break;
    }
  }
  return bytes_written;
}

IOWrapperSizeT ChunkedMpiByteWriteAtAll(IOWrapperFile fh, MPI_Comm comm,
                                        const char* buf,
                                        IOWrapperSizeT total_bytes,
                                        IOWrapperSizeT offset) {
  PreflightPositionedMpiRange(offset, total_bytes, "MPI_File_write_at_all");
  IOWrapperSizeT max_total_bytes =
      CollectiveMaxBytes(comm, total_bytes, "MPI_Allreduce write byte count");
  IOWrapperSizeT max_chunk =
      AgreedMaxMpiBytes(comm, "MPI_Allreduce write chunk limit");
  char dummy = '\0';
  IOWrapperSizeT bytes_written = 0;
  IOWrapperSizeT chunk_begin = 0;
  while (chunk_begin < max_total_bytes) {
    IOWrapperSizeT scheduled_bytes =
        std::min(max_chunk, max_total_bytes - chunk_begin);
    IOWrapperSizeT local_bytes = 0;
    const char* local_buf = &dummy;
    MPI_Offset mpi_offset = 0;
    if (chunk_begin < total_bytes) {
      local_bytes = std::min(scheduled_bytes, total_bytes - chunk_begin);
      local_buf = buf + CheckedSizeT(chunk_begin,
                                     "Write_any_type_at_all buffer offset");
      mpi_offset = CheckedMpiOffset(
          CheckedOffsetAdd(offset, chunk_begin, "MPI_File_write_at_all"),
          "MPI_File_write_at_all");
    }
    MPI_Status status;
    int mpi_error = MPI_File_write_at_all(
        fh, mpi_offset, const_cast<char*>(local_buf),
        static_cast<int>(local_bytes), MPI_BYTE, &status);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError("MPI_File_write_at_all", mpi_error);
    }
    bool valid = false;
    IOWrapperSizeT transferred = 0;
    if (mpi_error == MPI_SUCCESS) {
      transferred = TransferredBytes(
          &status, "MPI_File_write_at_all MPI_Get_count", &valid);
    }
    bool local_success =
        mpi_error == MPI_SUCCESS && valid && transferred == local_bytes;
    if (!CollectiveRoundSucceeded(comm, local_success,
                                  "MPI_Allreduce write result")) {
      return 0;
    }
    bytes_written = CheckedOffsetAdd(bytes_written, transferred,
                                     "MPI_File_write_at_all");
    chunk_begin = CheckedOffsetAdd(chunk_begin, scheduled_bytes,
                                   "MPI_File_write_at_all progress");
  }
  return bytes_written;
}

bool MpiFileIsMissing(int mpi_error) {
  int error_class = MPI_SUCCESS;
  int class_error = MPI_Error_class(mpi_error, &error_class);
  if (class_error != MPI_SUCCESS) {
    FatalMpiError("MPI_Error_class", class_error);
  }
  return error_class == MPI_ERR_NO_SUCH_FILE;
}

void DeleteExistingMpiFile(const char* fname, MPI_Comm comm) {
  int comm_rank = 0;
  int mpi_error = MPI_Comm_rank(comm, &comm_rank);
  if (mpi_error != MPI_SUCCESS) {
    FatalMpiError("MPI_Comm_rank", mpi_error);
  }
  if (comm_rank == 0) {
    mpi_error = MPI_File_delete(fname, MPI_INFO_NULL);
    if (mpi_error != MPI_SUCCESS && !MpiFileIsMissing(mpi_error)) {
      FatalMpiError("MPI_File_delete", mpi_error);
    }
  }
  mpi_error = MPI_Barrier(comm);
  if (mpi_error != MPI_SUCCESS) {
    FatalMpiError("MPI_Barrier before MPI_File_open", mpi_error);
  }
}
#endif

}  // namespace

namespace io_wrapper {

#if MPI_PARALLEL_ENABLED
void BroadcastBytes(void* buf, IOWrapperSizeT count, int root, MPI_Comm comm) {
  IOWrapperSizeT min_count = 0;
  IOWrapperSizeT max_count = 0;
  int mpi_error = MPI_Allreduce(&count, &min_count, 1, MPI_UINT64_T, MPI_MIN,
                                comm);
  if (mpi_error != MPI_SUCCESS) {
    FatalMpiError("MPI_Allreduce broadcast byte count", mpi_error);
  }
  mpi_error = MPI_Allreduce(&count, &max_count, 1, MPI_UINT64_T, MPI_MAX, comm);
  if (mpi_error != MPI_SUCCESS) {
    FatalMpiError("MPI_Allreduce broadcast byte count", mpi_error);
  }
  if (min_count != max_count) {
    FatalIOError("BroadcastBytes count must match on every communicator rank.");
  }

  char* byte_buf = reinterpret_cast<char*>(buf);
  IOWrapperSizeT max_chunk =
      AgreedMaxMpiBytes(comm, "MPI_Allreduce broadcast chunk limit");
  IOWrapperSizeT offset = 0;
  while (offset < count) {
    IOWrapperSizeT chunk_bytes = std::min(max_chunk, count - offset);
    mpi_error = MPI_Bcast(
        byte_buf + CheckedSizeT(offset, "BroadcastBytes buffer offset"),
        static_cast<int>(chunk_bytes), MPI_BYTE, root, comm);
    if (mpi_error != MPI_SUCCESS) {
      FatalMpiError("MPI_Bcast", mpi_error);
    }
    offset = CheckedOffsetAdd(offset, chunk_bytes, "MPI_Bcast progress");
  }
}
#endif

}  // namespace io_wrapper

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open(const char* fname, FileMode rw)
//! \brief wrapper for {MPI_File_open} versus {std::fopen} including error check
//! This function must not be called by multiple threads in shared memory parallel regions

int IOWrapper::Open(const char* fname, FileMode rw, bool use_serial_io) {
  const char* mode;
  switch (rw) {
    case FileMode::read:
      mode = "rb";
      break;
    case FileMode::write:
      mode = "wb";
      break;
    case FileMode::append:
      mode = "ab";
      break;
    default:
      return false;
  }

#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    int mpi_mode;
    switch (rw) {
      case FileMode::read:
        mpi_mode = MPI_MODE_RDONLY;
        break;
      case FileMode::write:
        mpi_mode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
        DeleteExistingMpiFile(fname, comm_);
        break;
      case FileMode::append:
        mpi_mode = MPI_MODE_WRONLY | MPI_MODE_APPEND;
        break;
      default:
        return false;
    }

    int mpi_error = MPI_File_open(comm_, fname, mpi_mode, MPI_INFO_NULL, &fh_);
    if (mpi_error != MPI_SUCCESS) {
      FatalMpiError("MPI_File_open", mpi_error);
    }
  } else {
    FILE* local_fh;
    if ((local_fh = std::fopen(fname, mode)) == nullptr) {
      perror("Error opening file");
      FatalIOError("File '" + std::string(fname) + "' could not be opened.");
    }
    fh_ = reinterpret_cast<IOWrapperFile>(local_fh);
  }
#else
  FILE* local_fh;
  if ((local_fh = std::fopen(fname, mode)) == nullptr) {
    perror("Error opening file");
    FatalIOError("File '" + std::string(fname) + "' could not be opened.");
  }
  fh_ = local_fh;
#endif

  return true;
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_read} versus {std::fread}.

std::size_t IOWrapper::Read_bytes(void* buf, IOWrapperSizeT size,
                                  IOWrapperSizeT cnt, bool use_serial_io) {
  if (size == 0 || cnt == 0) {
    return 0;
  }
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    IOWrapperSizeT total_bytes = CheckedByteCount(size, cnt, "Read_bytes");
    return CompletedElements(
        ChunkedMpiByteRead(fh_, reinterpret_cast<char*>(buf), total_bytes),
        size, "Read_bytes result");
  }
  return std::fread(buf, size, cnt, reinterpret_cast<FILE*>(fh_));
#else
  return std::fread(buf, size, cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_read_at} versus {std::fseek+std::fread}.

std::size_t IOWrapper::Read_bytes_at(void* buf, IOWrapperSizeT size,
                                     IOWrapperSizeT cnt, IOWrapperSizeT offset,
                                     bool use_serial_io) {
  if (size == 0 || cnt == 0) {
    return 0;
  }
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    IOWrapperSizeT total_bytes = CheckedByteCount(size, cnt, "Read_bytes_at");
    return CompletedElements(ChunkedMpiByteReadAt(
        fh_, reinterpret_cast<char*>(buf), total_bytes, offset),
        size, "Read_bytes_at result");
  }
  std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
  return std::fread(buf, size, cnt, reinterpret_cast<FILE*>(fh_));
#else
  std::fseek(fh_, offset, SEEK_SET);
  return std::fread(buf, size, cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_read_at_all} versus {std::fseek+std::fread}.

std::size_t IOWrapper::Read_bytes_at_all(void* buf, IOWrapperSizeT size,
                                         IOWrapperSizeT cnt,
                                         IOWrapperSizeT offset,
                                         bool use_serial_io) {
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    IOWrapperSizeT total_bytes =
        CheckedByteCount(size, cnt, "Read_bytes_at_all");
    char dummy = '\0';
    char* read_buf = total_bytes == 0 ? &dummy : reinterpret_cast<char*>(buf);
    IOWrapperSizeT read =
        ChunkedMpiByteReadAtAll(fh_, comm_, read_buf, total_bytes, offset);
    return size == 0 ? 0 : CompletedElements(read, size,
                                             "Read_bytes_at_all result");
  }
#endif
  if (size == 0 || cnt == 0) {
    return 0;
  }
#if MPI_PARALLEL_ENABLED
  std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
  return std::fread(buf, size, cnt, reinterpret_cast<FILE*>(fh_));
#else
  std::fseek(fh_, offset, SEEK_SET);
  return std::fread(buf, size, cnt, fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for reading Athena Reals.

std::size_t IOWrapper::Read_Reals(void* buf, IOWrapperSizeT cnt,
                                  bool use_serial_io) {
  return Read_bytes(buf, sizeof(Real), cnt, use_serial_io);
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for reading Athena Reals at an explicit offset.

std::size_t IOWrapper::Read_Reals_at(void* buf, IOWrapperSizeT cnt,
                                     IOWrapperSizeT offset,
                                     bool use_serial_io) {
  return Read_bytes_at(buf, sizeof(Real), cnt, offset, use_serial_io);
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for collectively reading Athena Reals at explicit offsets.

std::size_t IOWrapper::Read_Reals_at_all(void* buf, IOWrapperSizeT cnt,
                                         IOWrapperSizeT offset,
                                         bool use_serial_io) {
  return Read_bytes_at_all(buf, sizeof(Real), cnt, offset, use_serial_io);
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for writing any supported datatype.

std::size_t IOWrapper::Write_any_type(const void* buf, IOWrapperSizeT cnt,
                                      std::string datatype,
                                      bool use_serial_io) {
  std::size_t datasize = SizeOfDatatype(datatype);
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    IOWrapperSizeT total_bytes =
        CheckedByteCount(datasize, cnt, "Write_any_type");
    return CompletedElements(ChunkedMpiByteWrite(
        fh_, reinterpret_cast<const char*>(buf), total_bytes),
        datasize, "Write_any_type result");
  }
  std::size_t written =
      std::fwrite(buf, datasize, cnt, reinterpret_cast<FILE*>(fh_));
#else
  std::size_t written = std::fwrite(buf, datasize, cnt, fh_);
#endif
  if (written != cnt) {
    std::cerr << "Error writing data. Expected to write " << cnt
              << " elements, but wrote " << written << std::endl;
  }
  return written;
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for writing any supported datatype at an explicit offset.

std::size_t IOWrapper::Write_any_type_at(const void* buf, IOWrapperSizeT cnt,
                                         IOWrapperSizeT offset,
                                         std::string datatype,
                                         bool use_serial_io) {
  std::size_t datasize = SizeOfDatatype(datatype);
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    IOWrapperSizeT total_bytes =
        CheckedByteCount(datasize, cnt, "Write_any_type_at");
    return CompletedElements(ChunkedMpiByteWriteAt(
        fh_, reinterpret_cast<const char*>(buf), total_bytes, offset),
        datasize, "Write_any_type_at result");
  }
  std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
  std::size_t written =
      std::fwrite(buf, datasize, cnt, reinterpret_cast<FILE*>(fh_));
#else
  std::fseek(fh_, offset, SEEK_SET);
  std::size_t written = std::fwrite(buf, datasize, cnt, fh_);
#endif
  if (written != cnt) {
    std::cerr << "Error writing data. Expected to write " << cnt
              << " elements, but wrote " << written << std::endl;
  }
  return written;
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for collectively writing any supported datatype at explicit offsets.

std::size_t IOWrapper::Write_any_type_at_all(const void* buf,
                                             IOWrapperSizeT cnt,
                                             IOWrapperSizeT offset,
                                             std::string datatype,
                                             bool use_serial_io) {
  std::size_t datasize = SizeOfDatatype(datatype);
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    IOWrapperSizeT total_bytes =
        CheckedByteCount(datasize, cnt, "Write_any_type_at_all");
    char dummy = '\0';
    const char* write_buf =
        total_bytes == 0 ? &dummy : reinterpret_cast<const char*>(buf);
    return CompletedElements(ChunkedMpiByteWriteAtAll(
        fh_, comm_, write_buf, total_bytes, offset),
        datasize, "Write_any_type_at_all result");
  }
  std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
  std::size_t written =
      std::fwrite(buf, datasize, cnt, reinterpret_cast<FILE*>(fh_));
#else
  std::fseek(fh_, offset, SEEK_SET);
  std::size_t written = std::fwrite(buf, datasize, cnt, fh_);
#endif
  if (written != cnt) {
    std::cerr << "Error writing data. Expected to write " << cnt
              << " elements, but wrote " << written << std::endl;
  }
  return written;
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_close} versus {std::fclose}.

int IOWrapper::Close(bool use_serial_io) {
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    return MPI_File_close(&fh_);
  }
  return std::fclose(reinterpret_cast<FILE*>(fh_));
#else
  return std::fclose(fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_seek} versus {std::fseek}.

int IOWrapper::Seek(IOWrapperSizeT offset, bool use_serial_io) {
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    return MPI_File_seek(fh_, CheckedMpiOffset(offset, "MPI_File_seek"),
                         MPI_SEEK_SET);
  }
  return std::fseek(reinterpret_cast<FILE*>(fh_), offset, SEEK_SET);
#else
  return std::fseek(fh_, offset, SEEK_SET);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_get_position} versus {std::ftell}.

IOWrapperSizeT IOWrapper::GetPosition(bool use_serial_io) {
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    MPI_Offset position = 0;
    int mpi_error = MPI_File_get_position(fh_, &position);
    if (mpi_error != MPI_SUCCESS) {
      FatalMpiError("MPI_File_get_position", mpi_error);
    }
    if (position < 0) {
      FatalIOError("MPI_File_get_position returned a negative offset.");
    }
    return static_cast<IOWrapperSizeT>(position);
  }
  int64_t pos = ftell(reinterpret_cast<FILE*>(fh_));
  return pos;
#else
  int64_t pos = ftell(fh_);
  return pos;
#endif
}
