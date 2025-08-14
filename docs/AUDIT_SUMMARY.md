# AthenaK Documentation Audit Summary

## Date: 2024-08-13

## CRITICAL ISSUES FOUND AND FIXED

### 1. **MAJOR ERROR: Coordinate Systems (modules/coordinates.md)**
**Issue**: Documentation claimed AthenaK supports spherical polar and cylindrical coordinates
**Reality**: AthenaK ONLY supports:
- Cartesian coordinates (default)
- Kerr-Schild coordinates for black hole spacetimes (GR only)

**Status**: ✅ FIXED - Completely rewrote coordinates.md to accurately reflect only Cartesian and Kerr-Schild support

### 2. **Quickstart Input Files (quickstart.md)**
**Issue**: Referenced non-existent input files:
- `inputs/tests/shock_tube.athinput` - DOES NOT EXIST
- `inputs/tests/shock_tube_mhd.athinput` - DOES NOT EXIST
- `inputs/tests/blast.athinput` - DOES NOT EXIST

**Reality**: Correct files are:
- `inputs/hydro/sod.athinput` - Sod shock tube
- `inputs/mhd/bw.athinput` - Brio-Wu MHD shock
- `inputs/mhd/blast_mhd.athinput` - MHD blast wave

**Status**: ✅ FIXED - Updated all input file paths in quickstart.md

### 3. **Reconstruction Methods (hydro.md, mhd.md)**
**Issue**: Listed "ppm" as single option
**Reality**: AthenaK has:
- `ppm4` - 4th order PPM
- `ppmx` - Extremum-preserving PPM
- No generic "ppm"

**Status**: ✅ FIXED - Updated both hydro.md and mhd.md with correct reconstruction methods

### 4. **Ghost Zone Requirements**
**Issue**: Not documented that higher-order methods need more ghost zones
**Reality**: 
- PLM needs 3 ghost zones with FOFC
- PPM4/PPMX/WENOZ need 3 ghost zones (4 with FOFC)

**Status**: ✅ FIXED - Added ghost zone requirements table

### 5. **Mermaid Flowcharts**
**Issue**: 
- Flowcharts not rendering due to %%{init:} directives
- No click functionality for navigation
- Not working in dark mode

**Status**: ✅ FIXED - Removed problematic syntax, added click handlers, made theme-aware

## MODERATE ISSUES FIXED

### 6. **MHD Documentation**
- Fixed incorrect file references (e.g., "hlld_mhd.cpp" should be "hlld_mhd.hpp")
- Added note that Roe solver not implemented for MHD
- Added recent AMR div(B) preservation fixes
- Added FOFC documentation

### 7. **Hydro Documentation**
- Added viscosity and conductivity parameters
- Corrected reconstruction method names
- Added relativistic solver information

### 8. **Overview Page**
- All physics modules verified to exist
- Executable name confirmed as "athena"
- Fixed input file examples

## FILES AUDITED AND MODIFIED

### Completely Rewritten:
1. `modules/coordinates.md` - Major factual errors corrected

### Significantly Updated:
2. `modules/mhd.md` - Comprehensive update with accurate information
3. `modules/hydro.md` - Corrected reconstruction methods and parameters
4. `quickstart.md` - Fixed all input file paths

### Enhanced:
5. `flowcharts/runtime.md` - Added click functionality
6. `flowcharts/system_architecture.md` - Added click functionality  
7. `overview.md` - Added click functionality to flowchart
8. `_static/custom.js` - Made Mermaid diagrams theme-aware

## MODULES STILL TO AUDIT

The following documentation files have NOT been audited yet:
- modules/mesh.md
- modules/driver.md
- modules/tasklist.md
- modules/boundaries.md
- modules/outputs.md (partially checked, seems accurate)
- modules/reconstruction.md
- modules/riemann_solvers.md
- modules/eos.md
- modules/diffusion.md
- modules/srcterms.md
- modules/particles.md
- modules/radiation.md
- modules/z4c.md
- modules/dyn_grmhd.md
- modules/ion_neutral.md
- modules/shearing_box.md
- modules/pgen.md
- All example pages
- All reference documentation

## RECOMMENDATIONS

1. **URGENT**: Continue auditing remaining modules - there may be more inaccuracies
2. **Add CI Testing**: Create tests that verify documentation examples actually work
3. **Version Documentation**: Clearly state which version of AthenaK the docs apply to
4. **Add "Not Implemented" Section**: Document features from Athena++ not yet in AthenaK

## KEY FINDINGS

The documentation had significant inaccuracies, particularly around:
- Coordinate systems (completely wrong)
- Input file paths (non-existent files referenced)
- Some technical details (reconstruction methods, etc.)

This suggests the documentation may have been:
- Copied from Athena++ without verification
- Written speculatively without checking the actual code
- Not tested against the actual codebase

## NEXT STEPS

Continue systematic audit of remaining ~20 documentation files to ensure accuracy.