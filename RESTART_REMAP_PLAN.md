# Smart Restart File Remapping Implementation Plan
## Enabling Restart with Different MPI Rank Counts

### Executive Summary
This document outlines the implementation of smart restart file remapping in AthenaK, allowing simulations to restart with a different number of MPI ranks than the original run. This is essential for the turb_timed_amr problem where the domain is refined and restarted with 8x more ranks.

### Problem Statement
- **Current Issue**: When using `single_file_per_rank=true`, each rank expects its restart file in `rst/rank_XXXXXXXX/`
- **Use Case**: Run with N ranks → refine entire domain → restart with 8N ranks
- **Blocker**: New ranks cannot find restart files because rank IDs don't match

### Solution Architecture
Record the original rank layout in restart files and implement smart remapping that allows any rank to read data from any restart file based on MeshBlock ownership.

---

## Implementation Checklist

### Phase 1: Extend Restart File Format
Store original rank distribution information in restart headers.

#### [ ] Step 1.1: Update Header Writing (src/outputs/restart.cpp)
**Location**: `RestartOutput::WriteOutputFile()` around line 208-235

**Current Code**:
```cpp
// STEP 1: Root process writes header data
if (global_variable::my_rank == 0 || single_file_per_rank) {
    // output the input parameters
    resfile.Write_any_type(sbuf.c_str(), sbuf.size(), "byte", single_file_per_rank);
    // output Mesh information
    resfile.Write_any_type(&(pm->nmb_total), (sizeof(int)), "byte", single_file_per_rank);
    // ... other header data ...
}
```

**New Code** (add after existing header):
```cpp
// STEP 1b: Write rank distribution info for remapping
if (global_variable::my_rank == 0 || single_file_per_rank) {
    // Write original number of ranks
    resfile.Write_any_type(&(global_variable::nranks), sizeof(int), "byte", 
                           single_file_per_rank);
    // Write distribution of MeshBlocks per rank
    resfile.Write_any_type(pm->nmb_eachrank, global_variable::nranks*sizeof(int), 
                           "byte", single_file_per_rank);
    // Write starting global IDs per rank  
    resfile.Write_any_type(pm->gids_eachrank, global_variable::nranks*sizeof(int),
                           "byte", single_file_per_rank);
}
```

#### [ ] Step 1.2: Update Offset Calculations
**Location**: `RestartOutput::WriteOutputFile()` around line 298-310

**Current Calculation**:
```cpp
IOWrapperSizeT step1size = sbuf.size()*sizeof(char) + 3*sizeof(int) + 2*sizeof(Real) +
                          sizeof(RegionSize) + 2*sizeof(RegionIndcs);
```

**New Calculation**:
```cpp
IOWrapperSizeT step1size = sbuf.size()*sizeof(char) + 3*sizeof(int) + 2*sizeof(Real) +
                          sizeof(RegionSize) + 2*sizeof(RegionIndcs);
// Add size of new rank distribution data
IOWrapperSizeT step1bsize = sizeof(int) + 2*global_variable::nranks*sizeof(int);
IOWrapperSizeT offset_myrank = (step1size + step1bsize + step2size + step3size
                                + sizeof(IOWrapperSizeT));
```

---

### Phase 2: Detect and Map Restart Files

#### [ ] Step 2.1: Create Restart File Mapper Class
**New File**: `src/outputs/restart_mapper.hpp`

```cpp
#ifndef OUTPUTS_RESTART_MAPPER_HPP_
#define OUTPUTS_RESTART_MAPPER_HPP_

#include <string>
#include <map>
#include <vector>
#include "athena.hpp"
#include "outputs/io_wrapper.hpp"

class RestartFileMapper {
 public:
  RestartFileMapper() = default;
  ~RestartFileMapper();

  // Initialize mapping from available restart files
  void Initialize(const std::string& restart_path, bool single_file_per_rank);
  
  // Get source file info for a given MeshBlock
  struct SourceInfo {
    int orig_rank;           // Original rank that owned this MB
    std::string filepath;    // Path to restart file
    IOWrapperSizeT offset;   // Byte offset within file
  };
  SourceInfo GetSourceInfo(int mb_gid) const;
  
  // Read rank distribution from restart file
  void ReadRankDistribution(IOWrapper& resfile, bool single_file_per_rank);
  
  // Open/close file handles
  IOWrapper* OpenSourceFile(const std::string& filepath);
  void CloseAllFiles();

  // Original rank layout info
  int orig_nranks;
  int* orig_nmb_eachrank;
  int* orig_gids_eachrank;
  
 private:
  std::map<int, int> mb_to_orig_rank_;  // Maps MB GID to original rank
  std::map<std::string, IOWrapper*> open_files_;  // Cache of open file handles
  std::vector<std::string> available_files_;  // List of detected restart files
  
  void ScanAvailableFiles(const std::string& base_path);
  void BuildMBToRankMapping();
};

#endif // OUTPUTS_RESTART_MAPPER_HPP_
```

#### [ ] Step 2.2: Implement Restart File Mapper
**New File**: `src/outputs/restart_mapper.cpp`

```cpp
#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>
#include "outputs/restart_mapper.hpp"
#include "globals.hpp"

RestartFileMapper::~RestartFileMapper() {
  CloseAllFiles();
  delete[] orig_nmb_eachrank;
  delete[] orig_gids_eachrank;
}

void RestartFileMapper::Initialize(const std::string& restart_path, 
                                   bool single_file_per_rank) {
  if (single_file_per_rank) {
    ScanAvailableFiles(restart_path);
  } else {
    available_files_.push_back(restart_path);
  }
}

void RestartFileMapper::ScanAvailableFiles(const std::string& base_path) {
  // Extract base directory from path (before rank_XXXXXXXX)
  size_t rank_pos = base_path.find("/rank_");
  if (rank_pos == std::string::npos) {
    // Single file mode
    available_files_.push_back(base_path);
    return;
  }
  
  std::string base_dir = base_path.substr(0, rank_pos);
  size_t last_slash = base_path.rfind('/');
  std::string file_name = base_path.substr(last_slash + 1);
  
  // Scan for rank_* directories
  DIR* dir = opendir(base_dir.c_str());
  if (dir == nullptr) {
    std::cout << "### ERROR: Cannot open restart directory: " << base_dir << std::endl;
    return;
  }
  
  struct dirent* entry;
  while ((entry = readdir(dir)) != nullptr) {
    std::string dirname(entry->d_name);
    if (dirname.substr(0, 5) == "rank_") {
      std::string filepath = base_dir + "/" + dirname + "/" + file_name;
      struct stat st;
      if (stat(filepath.c_str(), &st) == 0) {
        available_files_.push_back(filepath);
      }
    }
  }
  closedir(dir);
  
  // Sort files by rank number for consistent ordering
  std::sort(available_files_.begin(), available_files_.end());
  
  if (global_variable::my_rank == 0) {
    std::cout << "Found " << available_files_.size() << " restart files" << std::endl;
  }
}

void RestartFileMapper::ReadRankDistribution(IOWrapper& resfile, 
                                            bool single_file_per_rank) {
  // Read original nranks
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read(&orig_nranks, 1, sizeof(int), "int", single_file_per_rank) != 1) {
      std::cout << "### FATAL ERROR: Failed to read original nranks from restart file\n";
      std::exit(EXIT_FAILURE);
    }
  }
  
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Bcast(&orig_nranks, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif

  // Allocate and read distribution arrays
  orig_nmb_eachrank = new int[orig_nranks];
  orig_gids_eachrank = new int[orig_nranks];
  
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read(orig_nmb_eachrank, orig_nranks, sizeof(int), "int", 
                     single_file_per_rank) != orig_nranks) {
      std::cout << "### FATAL ERROR: Failed to read nmb_eachrank from restart file\n";
      std::exit(EXIT_FAILURE);
    }
    if (resfile.Read(orig_gids_eachrank, orig_nranks, sizeof(int), "int",
                     single_file_per_rank) != orig_nranks) {
      std::cout << "### FATAL ERROR: Failed to read gids_eachrank from restart file\n";
      std::exit(EXIT_FAILURE);
    }
  }
  
#if MPI_PARALLEL_ENABLED
  if (!single_file_per_rank) {
    MPI_Bcast(orig_nmb_eachrank, orig_nranks, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(orig_gids_eachrank, orig_nranks, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif

  BuildMBToRankMapping();
}

void RestartFileMapper::BuildMBToRankMapping() {
  mb_to_orig_rank_.clear();
  for (int rank = 0; rank < orig_nranks; rank++) {
    int start_gid = orig_gids_eachrank[rank];
    int end_gid = start_gid + orig_nmb_eachrank[rank];
    for (int gid = start_gid; gid < end_gid; gid++) {
      mb_to_orig_rank_[gid] = rank;
    }
  }
}

RestartFileMapper::SourceInfo RestartFileMapper::GetSourceInfo(int mb_gid) const {
  SourceInfo info;
  
  // Find original rank that owned this MeshBlock
  auto it = mb_to_orig_rank_.find(mb_gid);
  if (it == mb_to_orig_rank_.end()) {
    std::cout << "### FATAL ERROR: MeshBlock " << mb_gid 
              << " not found in original rank mapping\n";
    std::exit(EXIT_FAILURE);
  }
  info.orig_rank = it->second;
  
  // Determine file path
  if (available_files_.size() == 1) {
    // Single file mode
    info.filepath = available_files_[0];
  } else {
    // Multi-file mode - map original rank to available file
    if (info.orig_rank < available_files_.size()) {
      info.filepath = available_files_[info.orig_rank];
    } else {
      std::cout << "### FATAL ERROR: No restart file for original rank " 
                << info.orig_rank << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  
  // Calculate offset will be done by caller based on data layout
  info.offset = 0;
  
  return info;
}

IOWrapper* RestartFileMapper::OpenSourceFile(const std::string& filepath) {
  auto it = open_files_.find(filepath);
  if (it != open_files_.end()) {
    return it->second;
  }
  
  IOWrapper* file = new IOWrapper();
  file->Open(filepath.c_str(), IOWrapper::FileMode::read, true);
  open_files_[filepath] = file;
  return file;
}

void RestartFileMapper::CloseAllFiles() {
  for (auto& pair : open_files_) {
    pair.second->Close(true);
    delete pair.second;
  }
  open_files_.clear();
}
```

---

### Phase 3: Update BuildTreeFromRestart

#### [ ] Step 3.1: Modify Header Reading
**File**: `src/mesh/build_tree.cpp`
**Location**: `BuildTreeFromRestart()` around line 319-376

**Add after line 326** (after calculating headersize):
```cpp
// Account for new rank distribution data in header
int temp_nranks;
if (single_file_per_rank) {
  // Read just nranks to calculate additional header size
  IOWrapperSizeT pos = resfile.GetPosition(single_file_per_rank);
  resfile.Seek(headeroffset + headersize, single_file_per_rank);
  resfile.Read(&temp_nranks, 1, sizeof(int), "int", single_file_per_rank);
  resfile.Seek(pos, single_file_per_rank);
} else if (global_variable::my_rank == 0) {
  IOWrapperSizeT pos = resfile.GetPosition(single_file_per_rank);
  resfile.Seek(headeroffset + headersize, single_file_per_rank);
  resfile.Read(&temp_nranks, 1, sizeof(int), "int", single_file_per_rank);
  resfile.Seek(pos, single_file_per_rank);
}

#if MPI_PARALLEL_ENABLED
if (!single_file_per_rank) {
  MPI_Bcast(&temp_nranks, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
#endif

IOWrapperSizeT rank_dist_size = sizeof(int) + 2*temp_nranks*sizeof(int);
```

#### [ ] Step 3.2: Pass Mapper to ProblemGenerator
**File**: `src/mesh/mesh.hpp`
**Add member variable** around line 200:
```cpp
class Mesh {
  // ... existing members ...
  RestartFileMapper* restart_mapper;  // For restart file remapping
  // ...
};
```

**File**: `src/mesh/build_tree.cpp`
**Initialize mapper** after reading header (around line 380):
```cpp
// Initialize restart file mapper if using single_file_per_rank
restart_mapper = nullptr;
if (single_file_per_rank) {
  restart_mapper = new RestartFileMapper();
  restart_mapper->Initialize(pin->GetString("job", "restart_file"), true);
  restart_mapper->ReadRankDistribution(resfile, true);
}
```

---

### Phase 4: Update ProblemGenerator Restart Constructor

#### [ ] Step 4.1: Modify Constructor Signature
**File**: `src/pgen/pgen.hpp`
**Update constructor** around line 70:
```cpp
class ProblemGenerator {
 public:
  ProblemGenerator(ParameterInput *pin, Mesh *pm);
  ProblemGenerator(ParameterInput *pin, Mesh *pm, IOWrapper resfile,
                  bool single_file_per_rank, RestartFileMapper* mapper = nullptr);
  // ...
};
```

#### [ ] Step 4.2: Implement Cross-File Reading
**File**: `src/pgen/pgen.cpp`
**Modify restart constructor** (lines 133-630):

**Replace lines 301-306** (offset calculation):
```cpp
// Old code:
int mygids = pm->gids_eachrank[global_variable::my_rank];
IOWrapperSizeT offset_myrank = headeroffset;
if (!single_file_per_rank) {
  offset_myrank += data_size_ * pm->gids_eachrank[global_variable::my_rank];
}
IOWrapperSizeT myoffset = offset_myrank;
```

**New code**:
```cpp
// Calculate offsets based on whether we're remapping
RestartFileMapper* mapper = pm->restart_mapper;
IOWrapperSizeT offset_myrank = headeroffset;
IOWrapperSizeT myoffset = offset_myrank;

if (mapper != nullptr && single_file_per_rank) {
  // Will calculate per-MB offsets during reading
  // Skip collective offset calculation
} else {
  // Traditional offset calculation
  int mygids = pm->gids_eachrank[global_variable::my_rank];
  if (!single_file_per_rank) {
    offset_myrank += data_size_ * mygids;
  }
  myoffset = offset_myrank;
}
```

**Replace MeshBlock reading loops** (e.g., lines 319-357 for hydro):
```cpp
if (phydro != nullptr) {
  Kokkos::realloc(ccin, nmb, nhydro, nout3, nout2, nout1);
  
  if (mapper != nullptr && single_file_per_rank) {
    // Remapped reading - each rank reads its MBs from appropriate files
    int my_gid_start = pm->gids_eachrank[global_variable::my_rank];
    
    for (int m = 0; m < pm->nmb_thisrank; m++) {
      int global_mb_id = my_gid_start + m;
      
      // Get source file info for this MeshBlock
      auto source = mapper->GetSourceInfo(global_mb_id);
      
      // Open source file if needed
      IOWrapper* source_file = mapper->OpenSourceFile(source.filepath);
      
      // Calculate offset for this MB in source file
      // Need to account for header, other ranks' data, etc.
      int mb_local_id = global_mb_id - mapper->orig_gids_eachrank[source.orig_rank];
      IOWrapperSizeT mb_offset = CalculateMBOffset(source.orig_rank, mb_local_id, 
                                                   data_size_, mapper);
      
      // Read this MeshBlock's data
      auto mbptr = Kokkos::subview(ccin, m, Kokkos::ALL, Kokkos::ALL, 
                                  Kokkos::ALL, Kokkos::ALL);
      int mbcnt = mbptr.size();
      
      if (source_file->Read_Reals_at(mbptr.data(), mbcnt, mb_offset, true) != mbcnt) {
        std::cout << "### FATAL ERROR: Failed to read MB " << global_mb_id 
                  << " from " << source.filepath << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  } else {
    // Original collective/non-collective reading code
    for (int m=0; m<noutmbs_max; ++m) {
      // ... existing code ...
    }
  }
  
  Kokkos::deep_copy(Kokkos::subview(phydro->u0, std::make_pair(0,nmb), Kokkos::ALL,
                    Kokkos::ALL, Kokkos::ALL, Kokkos::ALL), ccin);
}
```

#### [ ] Step 4.3: Add Helper Function for Offset Calculation
**Add to pgen.cpp**:
```cpp
namespace {
IOWrapperSizeT CalculateMBOffset(int orig_rank, int mb_local_id, 
                                 IOWrapperSizeT data_size_per_mb,
                                 RestartFileMapper* mapper) {
  // Calculate header size (must match restart.cpp writing)
  // This needs the actual parameter input size which varies
  // For now, assume it's stored or can be recalculated
  
  // Base offset includes:
  // 1. Input parameters
  // 2. Original header (mesh info)  
  // 3. Rank distribution info
  // 4. Logical locations and costs
  // 5. Step 3 data (z4c, tracker, turb)
  // 6. Data size marker
  
  // Then add offset for MeshBlocks before this one
  IOWrapperSizeT offset = /* calculated header offset */;
  offset += mb_local_id * data_size_per_mb;
  
  return offset;
}
}
```

---

### Phase 5: Update Main Startup Logic

#### [ ] Step 5.1: Fix File Path Detection
**File**: `src/main.cpp`
**Current logic** (lines 236-264) already handles rank directory remapping.
**Verify it works** with the new mapper system.

#### [ ] Step 5.2: Pass Mapper Through Construction
**File**: `src/main.cpp`
**Update ProblemGenerator construction** around line 333:
```cpp
if (pmesh->restart_mapper != nullptr) {
  pmesh->pgen = std::make_unique<ProblemGenerator>(pinput, pmesh, restartfile, 
                                                   single_file_per_rank,
                                                   pmesh->restart_mapper);
} else {
  pmesh->pgen = std::make_unique<ProblemGenerator>(pinput, pmesh, restartfile,
                                                   single_file_per_rank);
}
```

---

### Phase 6: Testing and Validation

#### [ ] Test 6.1: Single Rank Write/Read Test
1. Run with 1 rank, single_file_per_rank=true
2. Restart with 1 rank - verify identical results
3. Restart with 8 ranks - verify conservation

#### [ ] Test 6.2: Multi-Rank Scaling Test  
1. Run with 8 ranks, single_file_per_rank=true
2. Restart with 8 ranks - verify identical
3. Restart with 64 ranks - verify conservation

#### [ ] Test 6.3: Turbulence Driver Test
1. Run turb_timed_amr with 8 ranks to refinement trigger
2. Let AMR refine entire domain
3. Write restart with single_file_per_rank=true
4. Restart with 64 ranks
5. Verify:
   - All data loaded correctly
   - div(B) = 0 preserved
   - Turbulence driver state preserved
   - Simulation continues correctly

#### [ ] Test 6.4: Edge Cases
- [ ] Restart with fewer ranks than original (should fail gracefully)
- [ ] Missing restart files (clear error message)
- [ ] Corrupted restart file (detection and error)
- [ ] Mixed single/multi file modes

---

### Phase 7: Documentation and Cleanup

#### [ ] Step 7.1: Update User Documentation
- Add section to wiki about flexible restart capability
- Document single_file_per_rank best practices
- Add examples for common use cases

#### [ ] Step 7.2: Add Diagnostic Output
```cpp
if (global_variable::my_rank == 0 && mapper != nullptr) {
  std::cout << "=== Restart File Remapping ===" << std::endl;
  std::cout << "Original run: " << mapper->orig_nranks << " ranks" << std::endl;
  std::cout << "Current run: " << global_variable::nranks << " ranks" << std::endl;
  std::cout << "Found " << mapper->available_files_.size() << " restart files" << std::endl;
  std::cout << "Remapping MeshBlocks from old to new rank distribution..." << std::endl;
}
```

#### [ ] Step 7.3: Memory Cleanup
- Ensure RestartFileMapper destructor called in Mesh destructor
- Close all cached file handles
- Delete allocated arrays

---

## Implementation Order

### Priority 1 (Core Functionality)
1. Step 1.1-1.2: Extend restart file headers
2. Step 2.1-2.2: Create RestartFileMapper class
3. Step 3.1-3.2: Update BuildTreeFromRestart
4. Step 4.1-4.3: Update ProblemGenerator

### Priority 2 (Integration)
5. Step 5.1-5.2: Update main startup logic
6. Step 6.1-6.2: Basic testing

### Priority 3 (Robustness)
7. Step 6.3-6.4: Comprehensive testing
8. Step 7.1-7.3: Documentation and cleanup

---

## Risk Mitigation

### Backwards Compatibility
- Check for presence of rank distribution data in restart files
- Fall back to original behavior if not present
- Version flag in restart header for future changes

### Performance Considerations
- Cache open file handles to avoid repeated open/close
- Use MPI-IO collective operations where possible
- Consider memory-mapped files for large restarts

### Error Handling
- Clear error messages for common failure modes
- Graceful degradation when possible
- Diagnostic mode for debugging file mapping

---

## Success Criteria

1. **Functional**: Can restart turb_timed_amr with 8x more ranks
2. **Correct**: Conservation laws preserved, div(B)=0 maintained
3. **Performant**: No significant overhead vs. same-rank restart
4. **Robust**: Clear errors for edge cases, backwards compatible
5. **Documented**: Users can easily understand and use the feature