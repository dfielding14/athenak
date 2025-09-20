# Restart Remapping Implementation Checklist

## Quick Reference Implementation Steps

### 🔧 Core File Modifications

#### [ ] 1. src/outputs/restart.cpp
- [ ] Add rank distribution to header (after line 234)
  - Write `global_variable::nranks`
  - Write `nmb_eachrank` array  
  - Write `gids_eachrank` array
- [ ] Update offset calculations to account for new header size

#### [ ] 2. src/outputs/restart_mapper.hpp (NEW FILE)
- [ ] Create RestartFileMapper class definition
- [ ] Define SourceInfo struct
- [ ] Declare public methods: Initialize, GetSourceInfo, ReadRankDistribution
- [ ] Declare private members: file cache, MB-to-rank mapping

#### [ ] 3. src/outputs/restart_mapper.cpp (NEW FILE)
- [ ] Implement ScanAvailableFiles() - detect restart files in rank_* dirs
- [ ] Implement ReadRankDistribution() - read rank arrays from header
- [ ] Implement BuildMBToRankMapping() - map each MB to original rank
- [ ] Implement GetSourceInfo() - return file/offset for given MB
- [ ] Implement file handle caching

#### [ ] 4. src/mesh/build_tree.cpp
- [ ] Read rank distribution from restart header (line ~380)
- [ ] Create RestartFileMapper instance if single_file_per_rank
- [ ] Store mapper pointer in Mesh class
- [ ] Pass mapper to subsequent operations

#### [ ] 5. src/mesh/mesh.hpp
- [ ] Add `RestartFileMapper* restart_mapper` member
- [ ] Initialize to nullptr in constructor
- [ ] Delete in destructor

#### [ ] 6. src/pgen/pgen.hpp
- [ ] Add RestartFileMapper* parameter to restart constructor

#### [ ] 7. src/pgen/pgen.cpp
- [ ] Accept mapper in restart constructor
- [ ] For each physics module (hydro, mhd, etc):
  - [ ] Check if mapper exists and single_file_per_rank
  - [ ] If yes: Loop through local MBs, read from source files
  - [ ] If no: Use existing collective/non-collective read
- [ ] Implement CalculateMBOffset() helper function

#### [ ] 8. src/main.cpp
- [ ] Pass mapper to ProblemGenerator constructor (line ~333)

### ✅ Verification Steps

#### [ ] Compile and Link
- [ ] Add restart_mapper.cpp to CMakeLists.txt
- [ ] Include necessary headers
- [ ] Verify no undefined symbols

#### [ ] Unit Testing
- [ ] Test RestartFileMapper::ScanAvailableFiles with mock directory
- [ ] Test MB-to-rank mapping with known distribution
- [ ] Test file handle caching

#### [ ] Integration Testing  
- [ ] 1→1 rank restart (no remapping)
- [ ] 1→8 rank restart (expansion)
- [ ] 8→64 rank restart (8x expansion)
- [ ] Verify conservation (mass, momentum, energy)
- [ ] Verify div(B) = 0 preservation

#### [ ] Turbulence Driver Specific
- [ ] Preserve RNG state across restart
- [ ] Preserve forcing field
- [ ] Continue driver correctly after remap

### 📊 Test Matrix

| Original Ranks | Restart Ranks | Mode | Expected Result |
|---------------|--------------|------|-----------------|
| 1 | 1 | single_file | ✓ Works (no remap) |
| 1 | 8 | single_file | ✓ Remaps correctly |
| 8 | 8 | single_file | ✓ Works (no remap) |
| 8 | 64 | single_file | ✓ Remaps correctly |
| 64 | 8 | single_file | ✗ Error (fewer ranks) |
| 8 | 64 | multi_file | ✓ Uses original code |

### 🐛 Error Conditions to Handle

- [ ] Missing restart files → Clear error message with missing rank IDs
- [ ] Rank count mismatch → Inform user of original vs requested ranks
- [ ] Corrupted header → Detect size mismatches
- [ ] File I/O errors → Report which file failed
- [ ] Memory allocation failures → Clean shutdown

### 📝 Key Algorithms

```
MB_to_source_file(mb_gid):
  1. orig_rank = mb_to_orig_rank_map[mb_gid]
  2. file = available_files[orig_rank]  
  3. local_mb_id = mb_gid - orig_gids_eachrank[orig_rank]
  4. offset = header_size + local_mb_id * data_size_per_mb
  5. return (file, offset)
```

### 💡 Implementation Tips

1. **Start Small**: Test with tiny meshes first (2x2x2 blocks)
2. **Add Logging**: Verbose output during development
3. **Incremental Testing**: Test each physics module separately
4. **Use Debugger**: Set breakpoints in Read_Reals_at calls
5. **Verify Offsets**: Print calculated vs actual file positions

### 🚫 Common Pitfalls

- Forgetting to update ALL offset calculations
- Not broadcasting rank info to all MPI ranks
- File handle leaks from unclosed files
- Incorrect byte offset arithmetic
- Endianness issues between different systems

### 📈 Performance Optimizations (Future)

- [ ] Parallel file reading (multiple ranks read same file)
- [ ] Memory-mapped I/O for large files
- [ ] Prefetch next MB while processing current
- [ ] MPI collective I/O where beneficial

### ✨ Success Indicators

1. `div(B) = 0` to machine precision after restart
2. Identical evolution with/without restart at same rank count
3. No memory leaks (valgrind clean)
4. Performance within 10% of same-rank restart
5. Clear user messages for all error cases