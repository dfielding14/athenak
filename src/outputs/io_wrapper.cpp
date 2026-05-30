//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file io_wrapper.cpp
//! \brief functions that provide wrapper for MPI-IO versus serial input/output

#include "io_wrapper.hpp"

#include <sys/types.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "mpi_utils.hpp"

namespace {

constexpr const char* kTestMaxMpiBytesEnv = "ATHENAK_TEST_MAX_MPI_BYTES";

[[noreturn]] void FatalIOError(const std::string& message) {
  mpi_utils::AbortWorld(std::string("### FATAL ERROR in ") + __FILE__ +
                        " at line " + std::to_string(__LINE__) + "\n" +
                        message);
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

off_t CheckedSerialOffset(IOWrapperSizeT offset, const char* context) {
  if constexpr (std::numeric_limits<off_t>::digits <
                std::numeric_limits<IOWrapperSizeT>::digits) {
    if (offset > static_cast<IOWrapperSizeT>(std::numeric_limits<off_t>::max())) {
      FatalIOError(std::string(context) + " exceeds off_t range.");
    }
  }
  return static_cast<off_t>(offset);
}

int SerialSeek(FILE* file, IOWrapperSizeT offset, const char* context) {
  return ::fseeko(file, CheckedSerialOffset(offset, context), SEEK_SET);
}

void RequireSerialSeek(FILE* file, IOWrapperSizeT offset, const char* context) {
  if (SerialSeek(file, offset, context) != 0) {
    FatalIOError(std::string(context) + " serial seek failed.");
  }
}

void PreflightPositionedSerialRange(IOWrapperSizeT offset,
                                    IOWrapperSizeT total_bytes,
                                    const char* context) {
  if (total_bytes == 0) {
    return;
  }
  IOWrapperSizeT end_offset = CheckedOffsetAdd(offset, total_bytes - 1, context);
  CheckedSerialOffset(end_offset, context);
}

IOWrapperSizeT SerialPosition(FILE* file, const char* context) {
  off_t position = ::ftello(file);
  if (position < 0) {
    FatalIOError(std::string(context) + " returned a negative offset.");
  }
  return static_cast<IOWrapperSizeT>(position);
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
std::string MpiIOContext(const char* operation, const std::string &path) {
  return std::string(operation) + " for IOWrapper path '" + path + "'";
}

std::string MpiErrorMessage(const std::string &context, int mpi_error) {
  return context + " failed with MPI error: " +
         mpi_utils::MpiErrorString(mpi_error);
}

void PrintMpiError(const std::string &context, int mpi_error) {
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

IOWrapperSizeT TransferredBytes(MPI_Status* status, const std::string &context,
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
                              const std::string &context) {
  int local_failed = local_success ? 0 : 1;
  int any_failed = 0;
  int mpi_error = MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX,
                                comm);
  mpi_utils::CheckMpi(mpi_error, context.c_str());
  return any_failed == 0;
}

IOWrapperSizeT CollectiveMaxBytes(MPI_Comm comm, IOWrapperSizeT local_bytes,
                                  const std::string &context) {
  IOWrapperSizeT max_bytes = 0;
  int mpi_error = MPI_Allreduce(&local_bytes, &max_bytes, 1, MPI_UINT64_T,
                                MPI_MAX, comm);
  mpi_utils::CheckMpi(mpi_error, context.c_str());
  return max_bytes;
}

IOWrapperSizeT AgreedMaxMpiBytes(MPI_Comm comm, const std::string &context) {
  IOWrapperSizeT local_bytes = GetConfiguredMaxMpiBytes();
  IOWrapperSizeT min_bytes = 0;
  IOWrapperSizeT max_bytes = 0;
  int mpi_error = MPI_Allreduce(&local_bytes, &min_bytes, 1, MPI_UINT64_T,
                                MPI_MIN, comm);
  mpi_utils::CheckMpi(mpi_error, context.c_str());
  mpi_error = MPI_Allreduce(&local_bytes, &max_bytes, 1, MPI_UINT64_T, MPI_MAX,
                            comm);
  mpi_utils::CheckMpi(mpi_error, context.c_str());
  if (min_bytes != max_bytes) {
    FatalIOError(std::string(kTestMaxMpiBytesEnv) +
                 " must have the same value on every communicator rank.");
  }
  return min_bytes;
}

IOWrapperSizeT ChunkedMpiByteRead(IOWrapperFile fh, char* buf,
                                  IOWrapperSizeT total_bytes,
                                  const std::string &path) {
  const std::string operation = MpiIOContext("MPI_File_read", path);
  IOWrapperSizeT bytes_read = 0;
  IOWrapperSizeT max_chunk = GetConfiguredMaxMpiBytes();
  while (bytes_read < total_bytes) {
    IOWrapperSizeT chunk_bytes = std::min(max_chunk, total_bytes - bytes_read);
    MPI_Status status;
    int mpi_error = MPI_File_read(
        fh, buf + CheckedSizeT(bytes_read, "Read_bytes buffer offset"),
        static_cast<int>(chunk_bytes), MPI_BYTE, &status);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError(operation, mpi_error);
      return 0;
    }
    bool valid = false;
    IOWrapperSizeT transferred = TransferredBytes(
        &status, MpiIOContext("MPI_File_read MPI_Get_count", path), &valid);
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
                                    IOWrapperSizeT offset,
                                    const std::string &path) {
  const std::string operation = MpiIOContext("MPI_File_read_at", path);
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
      PrintMpiError(operation, mpi_error);
      return 0;
    }
    bool valid = false;
    IOWrapperSizeT transferred = TransferredBytes(
        &status, MpiIOContext("MPI_File_read_at MPI_Get_count", path), &valid);
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
                                       IOWrapperSizeT offset,
                                       const std::string &path) {
  const std::string operation = MpiIOContext("MPI_File_read_at_all", path);
  PreflightPositionedMpiRange(offset, total_bytes, "MPI_File_read_at_all");
  IOWrapperSizeT max_total_bytes =
      CollectiveMaxBytes(comm, total_bytes,
                         MpiIOContext("MPI_Allreduce read byte count", path));
  IOWrapperSizeT max_chunk =
      AgreedMaxMpiBytes(comm, MpiIOContext("MPI_Allreduce read chunk limit", path));
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
      PrintMpiError(operation, mpi_error);
    }
    bool valid = false;
    IOWrapperSizeT transferred = 0;
    if (mpi_error == MPI_SUCCESS) {
      transferred = TransferredBytes(
          &status, MpiIOContext("MPI_File_read_at_all MPI_Get_count", path), &valid);
    }
    bool local_success =
        mpi_error == MPI_SUCCESS && valid && transferred == local_bytes;
    if (!CollectiveRoundSucceeded(comm, local_success,
                                  MpiIOContext("MPI_Allreduce read result", path))) {
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
                                   IOWrapperSizeT total_bytes,
                                   const std::string &path) {
  const std::string operation = MpiIOContext("MPI_File_write", path);
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
      PrintMpiError(operation, mpi_error);
      return 0;
    }
    bool valid = false;
    IOWrapperSizeT transferred = TransferredBytes(
        &status, MpiIOContext("MPI_File_write MPI_Get_count", path), &valid);
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
                                     IOWrapperSizeT offset,
                                     const std::string &path) {
  const std::string operation = MpiIOContext("MPI_File_write_at", path);
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
      PrintMpiError(operation, mpi_error);
      return 0;
    }
    bool valid = false;
    IOWrapperSizeT transferred = TransferredBytes(
        &status, MpiIOContext("MPI_File_write_at MPI_Get_count", path), &valid);
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
                                        IOWrapperSizeT offset,
                                        const std::string &path) {
  const std::string operation = MpiIOContext("MPI_File_write_at_all", path);
  PreflightPositionedMpiRange(offset, total_bytes, "MPI_File_write_at_all");
  IOWrapperSizeT max_total_bytes =
      CollectiveMaxBytes(comm, total_bytes,
                         MpiIOContext("MPI_Allreduce write byte count", path));
  IOWrapperSizeT max_chunk =
      AgreedMaxMpiBytes(comm, MpiIOContext("MPI_Allreduce write chunk limit", path));
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
      PrintMpiError(operation, mpi_error);
    }
    bool valid = false;
    IOWrapperSizeT transferred = 0;
    if (mpi_error == MPI_SUCCESS) {
      transferred = TransferredBytes(
          &status, MpiIOContext("MPI_File_write_at_all MPI_Get_count", path), &valid);
    }
    bool local_success =
        mpi_error == MPI_SUCCESS && valid && transferred == local_bytes;
    if (!CollectiveRoundSucceeded(comm, local_success,
                                  MpiIOContext("MPI_Allreduce write result", path))) {
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
  mpi_utils::CheckMpi(class_error, "MPI_Error_class");
  return error_class == MPI_ERR_NO_SUCH_FILE;
}

void DeleteExistingMpiFile(const char* fname, MPI_Comm comm) {
  int comm_rank = 0;
  int mpi_error = MPI_Comm_rank(comm, &comm_rank);
  mpi_utils::CheckMpi(mpi_error,
                      "MPI_Comm_rank for IOWrapper publication communicator");
  if (comm_rank == 0) {
    mpi_error = MPI_File_delete(fname, MPI_INFO_NULL);
    if (mpi_error != MPI_SUCCESS && !MpiFileIsMissing(mpi_error)) {
      mpi_utils::CheckMpi(mpi_error,
                          "MPI_File_delete for IOWrapper replacement");
    }
  }
  mpi_error = MPI_Barrier(comm);
  mpi_utils::CheckMpi(
      mpi_error, "MPI_Barrier before MPI_File_open for IOWrapper publication");
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
  mpi_utils::CheckMpi(mpi_error, "MPI_Allreduce broadcast byte count");
  mpi_error = MPI_Allreduce(&count, &max_count, 1, MPI_UINT64_T, MPI_MAX, comm);
  mpi_utils::CheckMpi(mpi_error, "MPI_Allreduce broadcast byte count");
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
    mpi_utils::CheckMpi(mpi_error, "MPI_Bcast for IOWrapper byte broadcast");
    offset = CheckedOffsetAdd(offset, chunk_bytes,
                              "MPI_Bcast for IOWrapper byte-broadcast progress");
  }
}
#endif

}  // namespace io_wrapper

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open(const char* fname, FileMode rw)
//! \brief wrapper for {MPI_File_open} versus {std::fopen} including error check
//! This function must not be called by multiple threads in shared memory parallel regions

int IOWrapper::Open(const char* fname, FileMode rw, bool use_serial_io) {
  path_ = fname;
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
    std::string open_context =
        "MPI_File_open for IOWrapper path '" + std::string(fname) + "'";
    mpi_utils::CheckMpi(mpi_error, open_context.c_str());
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
        ChunkedMpiByteRead(fh_, reinterpret_cast<char*>(buf), total_bytes, path_),
        size, "Read_bytes result");
  }
  return std::fread(buf, CheckedSizeT(size, "Read_bytes element size"),
                    CheckedSizeT(cnt, "Read_bytes element count"),
                    reinterpret_cast<FILE*>(fh_));
#else
  return std::fread(buf, CheckedSizeT(size, "Read_bytes element size"),
                    CheckedSizeT(cnt, "Read_bytes element count"), fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_read_at} versus checked {fseeko+std::fread}.

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
        fh_, reinterpret_cast<char*>(buf), total_bytes, offset, path_),
        size, "Read_bytes_at result");
  }
  IOWrapperSizeT total_bytes = CheckedByteCount(size, cnt, "Read_bytes_at");
  PreflightPositionedSerialRange(offset, total_bytes, "Read_bytes_at");
  RequireSerialSeek(reinterpret_cast<FILE*>(fh_), offset, "Read_bytes_at");
  return std::fread(buf, CheckedSizeT(size, "Read_bytes_at element size"),
                    CheckedSizeT(cnt, "Read_bytes_at element count"),
                    reinterpret_cast<FILE*>(fh_));
#else
  IOWrapperSizeT total_bytes = CheckedByteCount(size, cnt, "Read_bytes_at");
  PreflightPositionedSerialRange(offset, total_bytes, "Read_bytes_at");
  RequireSerialSeek(fh_, offset, "Read_bytes_at");
  return std::fread(buf, CheckedSizeT(size, "Read_bytes_at element size"),
                    CheckedSizeT(cnt, "Read_bytes_at element count"), fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_read_at_all} versus checked {fseeko+std::fread}.

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
        ChunkedMpiByteReadAtAll(fh_, comm_, read_buf, total_bytes, offset, path_);
    return size == 0 ? 0 : CompletedElements(read, size,
                                             "Read_bytes_at_all result");
  }
#endif
  if (size == 0 || cnt == 0) {
    return 0;
  }
#if MPI_PARALLEL_ENABLED
  IOWrapperSizeT total_bytes = CheckedByteCount(size, cnt, "Read_bytes_at_all");
  PreflightPositionedSerialRange(offset, total_bytes, "Read_bytes_at_all");
  RequireSerialSeek(reinterpret_cast<FILE*>(fh_), offset, "Read_bytes_at_all");
  return std::fread(buf, CheckedSizeT(size, "Read_bytes_at_all element size"),
                    CheckedSizeT(cnt, "Read_bytes_at_all element count"),
                    reinterpret_cast<FILE*>(fh_));
#else
  IOWrapperSizeT total_bytes = CheckedByteCount(size, cnt, "Read_bytes_at_all");
  PreflightPositionedSerialRange(offset, total_bytes, "Read_bytes_at_all");
  RequireSerialSeek(fh_, offset, "Read_bytes_at_all");
  return std::fread(buf, CheckedSizeT(size, "Read_bytes_at_all element size"),
                    CheckedSizeT(cnt, "Read_bytes_at_all element count"), fh_);
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
        fh_, reinterpret_cast<const char*>(buf), total_bytes, path_),
        datasize, "Write_any_type result");
  }
  std::size_t written = std::fwrite(
      buf, datasize, CheckedSizeT(cnt, "Write_any_type element count"),
      reinterpret_cast<FILE*>(fh_));
#else
  std::size_t written =
      std::fwrite(buf, datasize, CheckedSizeT(cnt, "Write_any_type element count"), fh_);
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
        fh_, reinterpret_cast<const char*>(buf), total_bytes, offset, path_),
        datasize, "Write_any_type_at result");
  }
  IOWrapperSizeT total_bytes = CheckedByteCount(datasize, cnt, "Write_any_type_at");
  PreflightPositionedSerialRange(offset, total_bytes, "Write_any_type_at");
  RequireSerialSeek(reinterpret_cast<FILE*>(fh_), offset, "Write_any_type_at");
  std::size_t written = std::fwrite(
      buf, datasize, CheckedSizeT(cnt, "Write_any_type_at element count"),
      reinterpret_cast<FILE*>(fh_));
#else
  IOWrapperSizeT total_bytes = CheckedByteCount(datasize, cnt, "Write_any_type_at");
  PreflightPositionedSerialRange(offset, total_bytes, "Write_any_type_at");
  RequireSerialSeek(fh_, offset, "Write_any_type_at");
  std::size_t written = std::fwrite(
      buf, datasize, CheckedSizeT(cnt, "Write_any_type_at element count"), fh_);
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
        fh_, comm_, write_buf, total_bytes, offset, path_),
        datasize, "Write_any_type_at_all result");
  }
  IOWrapperSizeT total_bytes =
      CheckedByteCount(datasize, cnt, "Write_any_type_at_all");
  PreflightPositionedSerialRange(offset, total_bytes, "Write_any_type_at_all");
  RequireSerialSeek(reinterpret_cast<FILE*>(fh_), offset, "Write_any_type_at_all");
  std::size_t written = std::fwrite(
      buf, datasize, CheckedSizeT(cnt, "Write_any_type_at_all element count"),
      reinterpret_cast<FILE*>(fh_));
#else
  IOWrapperSizeT total_bytes =
      CheckedByteCount(datasize, cnt, "Write_any_type_at_all");
  PreflightPositionedSerialRange(offset, total_bytes, "Write_any_type_at_all");
  RequireSerialSeek(fh_, offset, "Write_any_type_at_all");
  std::size_t written = std::fwrite(
      buf, datasize, CheckedSizeT(cnt, "Write_any_type_at_all element count"), fh_);
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
    int mpi_error = MPI_File_close(&fh_);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError(MpiIOContext("MPI_File_close", path_), mpi_error);
    }
    return mpi_error;
  }
  return std::fclose(reinterpret_cast<FILE*>(fh_));
#else
  return std::fclose(fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_seek} versus checked {fseeko}.

int IOWrapper::Seek(IOWrapperSizeT offset, bool use_serial_io) {
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    int mpi_error = MPI_File_seek(fh_, CheckedMpiOffset(offset, "MPI_File_seek"),
                                  MPI_SEEK_SET);
    if (mpi_error != MPI_SUCCESS) {
      PrintMpiError(MpiIOContext("MPI_File_seek", path_), mpi_error);
    }
    return mpi_error;
  }
  return SerialSeek(reinterpret_cast<FILE*>(fh_), offset, "Seek");
#else
  return SerialSeek(fh_, offset, "Seek");
#endif
}

//----------------------------------------------------------------------------------------
//! \brief wrapper for {MPI_File_get_position} versus checked {ftello}.

IOWrapperSizeT IOWrapper::GetPosition(bool use_serial_io) {
#if MPI_PARALLEL_ENABLED
  if (!use_serial_io) {
    MPI_Offset position = 0;
    int mpi_error = MPI_File_get_position(fh_, &position);
    std::string context = MpiIOContext("MPI_File_get_position", path_);
    mpi_utils::CheckMpi(mpi_error, context.c_str());
    if (position < 0) {
      FatalIOError(context + " returned a negative offset.");
    }
    return static_cast<IOWrapperSizeT>(position);
  }
  return SerialPosition(reinterpret_cast<FILE*>(fh_), "GetPosition");
#else
  return SerialPosition(fh_, "GetPosition");
#endif
}
