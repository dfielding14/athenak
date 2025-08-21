# AthenaK Documentation Complete Overhaul Plan

## Current State Assessment

### Critical Errors Found
1. **Factual Errors**: Documentation claims HDF5 support - **DOES NOT EXIST**
2. **Missing Content**: Only 2 modules documented out of 20+ actual modules
3. **Broken Graphics**: Mermaid flowcharts rendering as raw text
4. **Incomplete Coverage**: Missing documentation for reconstruction, Riemann solvers, source terms, etc.
5. **Wrong Defaults**: Failed to document actual default output format (none - must be specified)
6. **Missing Parameters**: `single_file_per_rank` and other critical flags not documented

### Actual Codebase Structure (Verified)
```
src/
├── Core Infrastructure
│   ├── main.cpp, athena.hpp, globals.hpp
│   ├── driver/ - Time evolution control
│   ├── mesh/ - AMR mesh infrastructure
│   └── tasklist/ - Task-based execution
├── Physics Modules (7 total)
│   ├── hydro/ - Hydrodynamics
│   ├── mhd/ - Magnetohydrodynamics
│   ├── radiation/ - Radiation transport
│   ├── z4c/ - Numerical relativity
│   ├── dyn_grmhd/ - GR MHD in dynamical spacetimes
│   ├── ion-neutral/ - Two-fluid MHD
│   └── particles/ - Lagrangian particles
├── Numerical Methods
│   ├── reconstruct/ - PLM, PPM, WENOZ, DC
│   ├── eos/ - Ideal, isothermal, primitive solvers
│   └── diffusion/ - Viscosity, resistivity, conduction
├── Support Systems
│   ├── bvals/ - Boundary values and MPI
│   ├── coordinates/ - Cartesian, spherical, ADM
│   ├── outputs/ - 14 output formats (NO HDF5!)
│   ├── srcterms/ - Source terms, turbulence driving
│   └── shearing_box/ - Orbital advection
└── Problem Generators
    └── pgen/ - 50+ test problems and setups
```

## Documentation Overhaul Plan

### Phase 1: Fix Critical Infrastructure (Day 1)

#### 1.1 Fix Mermaid Flowcharts (ESSENTIAL - MUST WORK)
**Problem**: Flowcharts displaying as raw text, not rendering as diagrams
**Root Cause**: Sphinx/MyST not properly configured for Mermaid

**Complete Fix Process**:
```bash
# Step 1: Install proper dependencies
pip install sphinxcontrib-mermaid
npm install -g @mermaid-js/mermaid-cli

# Step 2: Update conf.py
extensions = [
    'sphinxcontrib.mermaid',
    # ... other extensions
]

# Mermaid configuration
mermaid_version = "10.6.0"
mermaid_init_js = "mermaid.initialize({startOnLoad:true, theme:'dark'});"

# Step 3: Fix flowchart syntax in markdown files
```

**Correct Mermaid Syntax**:
````markdown
```{mermaid}
flowchart TD
    Start([Start]) --> Init[Initialize MPI]
    Init --> Kokkos[Initialize Kokkos]
    Kokkos --> Driver[Create Driver]
    Driver --> Loop{Time Loop}
    Loop -->|t < tlim| Execute[Execute Tasks]
    Execute --> Loop
    Loop -->|t >= tlim| End([End])
```
````

**Alternative: Pre-render to SVG**:
```bash
# For each mermaid diagram:
mmdc -i diagram.mmd -o diagram.svg -t dark -b transparent
# Then embed SVG in documentation
```

**Testing Process**:
1. Build docs: `make html`
2. Open in browser
3. Verify ALL flowcharts render
4. Check clickable links work
5. Test in Chrome, Firefox, Safari

#### 1.2 Create Master Flowcharts (ESSENTIAL)
- [ ] **System Architecture Flowchart**: Shows all modules and data flow
- [ ] **Runtime Execution Flowchart**: Accurate main.cpp → driver flow
- [ ] **Task Dependency Graph**: How tasks execute in order
- [ ] **AMR Refinement Flowchart**: Mesh refinement process
- [ ] **I/O Pipeline Flowchart**: Output system flow

### Phase 2: Systematic Content Verification (Day 2-3)

#### 2.1 Verification Checklist Template
For EVERY documentation page, verify:
```markdown
## Verification Checklist for [Page Name]

### Source Code Verification
- [ ] Open corresponding source files
- [ ] Read actual implementation
- [ ] List all functions/classes
- [ ] Check default values
- [ ] Verify parameter names

### Cross-Reference Check
- [ ] All mentioned files exist
- [ ] All parameter names correct
- [ ] All class/function names accurate
- [ ] Default values match code
- [ ] Available options complete

### Testing
- [ ] Build example with documented parameters
- [ ] Run and verify behavior
- [ ] Check output matches description
```

#### 2.2 Module Documentation Template
```markdown
# Module: [Name]

## Source Location
`src/[directory]/`

## Purpose
[1-2 sentences from actual code comments]

## Files
| File | Purpose | Verified |
|------|---------|----------|
| [file.cpp] | [from code] | ✓/✗ |

## Configuration Parameters
| Parameter | Type | Default | Source File | Line |
|-----------|------|---------|-------------|------|
| [name] | [type] | [actual default] | [file] | [L###] |

## Dependencies
[List from actual includes]

## Verification Status
- Source reviewed: [date]
- Tested: [date]
- Reviewer: [name]
```

### Phase 2.5: Complete Input Parameter Documentation (CRITICAL)

#### Input Block Documentation Template
```markdown
# Input Block: <[block_name]>

## Module/File
- Source file: `src/[path]/[file].cpp`
- Line numbers: L[start]-L[end]
- Module: [which module uses this]

## Parameters
| Parameter | Type | Default | Required | Description | Source Ref |
|-----------|------|---------|----------|-------------|------------|
| [name] | int/Real/string/bool | [value] | Y/N | [what it does] | [file:line] |

## Example Usage
```ini
<[block_name]>
param1 = value1  # comment
param2 = value2
```

## Validation Rules
- [Any constraints or valid ranges]
- [Dependencies on other parameters]

## Related Blocks
- Links to other blocks that interact
```

#### Complete List of Input Blocks (from codebase scan)
```
Essential Blocks:
<job> - Job control parameters
<time> - Time integration control  
<mesh> - Mesh configuration
<meshblock> - MeshBlock sizing
<hydro> - Hydrodynamics parameters
<mhd> - MHD parameters
<eos> - Equation of state
<problem> - Problem-specific parameters
<coord> - Coordinate system
<refinement> - AMR control

Output Blocks:
<output1> through <output99> - Output configuration

Physics Blocks:
<radiation> - Radiation parameters
<z4c> - Numerical relativity
<particles> - Particle parameters
<ion-neutral> - Two-fluid parameters
<dyn_grmhd> - GR MHD parameters

Source/Diffusion Blocks:
<srcterms> - Source terms
<turb> - Turbulence driving
<viscosity> - Viscous diffusion
<resistivity> - Ohmic diffusion
<conduction> - Thermal conduction

Special Blocks:
<shearing_box> - Shearing box parameters
<units> - Unit system
<restart> - Restart configuration
```

### Phase 3: Complete Module Documentation (Day 4-5)

#### 3.1 Required Module Pages (20 total)
```
Core (4):
- [ ] mesh.md - AMR infrastructure
- [ ] driver.md - Time integration
- [ ] tasklist.md - Task execution
- [ ] coordinates.md - Coordinate systems

Physics (7):
- [ ] hydro.md - Hydrodynamics
- [ ] mhd.md - Magnetohydrodynamics  
- [ ] radiation.md - Radiation transport
- [ ] z4c.md - Numerical relativity
- [ ] dyn_grmhd.md - GR MHD
- [ ] ion_neutral.md - Two-fluid
- [ ] particles.md - Lagrangian particles

Numerical (4):
- [ ] reconstruction.md - PLM, PPM, WENOZ, DC
- [ ] riemann_solvers.md - All 20+ solvers
- [ ] eos.md - Equations of state
- [ ] diffusion.md - Physical diffusion

Support (5):
- [ ] outputs.md - All 14 formats (NO HDF5!)
- [ ] boundaries.md - Boundary conditions
- [ ] srcterms.md - Source terms
- [ ] shearing_box.md - Orbital advection
- [ ] pgen.md - Problem generators
```

#### 3.2 Output Format Documentation (CRITICAL FIX)
```markdown
## Actual Output Formats (from outputs.cpp)

### Mesh Data Outputs
- `vtk` - VTK legacy format (ParaView/VisIt compatible)
- `bin` - Binary format (native, fast)
- `rst` - Restart files (exact state preservation)
- `tab` - Formatted ASCII tables
- `cbin` - Coarsened binary (reduced resolution)

### Particle Outputs (if particles enabled)
- `pvtk` - Particle VTK format
- `trk` - Tracked particle trajectories
- `ppd` - Particle position dumps
- `prst` - Particle restart files

### Diagnostic Outputs
- `hst` - History files (time series)
- `log` - Event counters
- `pdf` - Probability density functions

### Key Parameters
- `single_file_per_rank` - Each MPI rank writes separately
- `ghost_zones` - Include ghost cells
- `coarsen_factor` - For cbin output

### NO HDF5 SUPPORT - Must remove all references
```

### Phase 4: File-by-File Reference (Day 6)

#### 4.1 Complete File Inventory
For each directory in src/:
```python
Script to generate:
1. List all .cpp and .hpp files
2. Extract first comment block
3. List public functions/classes
4. Generate markdown table
5. Verify against existing docs
```

#### 4.2 File Documentation Format
```markdown
## Directory: src/[name]/

### Purpose
[From directory-level README or main header]

### Files
| File | Lines | Purpose | Classes/Functions |
|------|-------|---------|-------------------|
| [name.cpp] | [###] | [extracted] | [list] |

### Dependencies
[Generated from includes]
```

### Phase 5: Testing and Validation (Day 7)

#### 5.1 Documentation Testing Protocol
```bash
#!/bin/bash
# For each documented example:

# 1. Copy example code/config
cp docs/examples/[example] test/

# 2. Build with documented parameters
cmake -B build [documented flags]
make -j8

# 3. Run with documented input
./athena -i [documented.athinput]

# 4. Verify output format
ls *.vtk *.bin  # Check actual output

# 5. Document discrepancies
echo "Expected: [doc], Actual: [result]" >> verification.log
```

#### 5.2 Automated Verification
```python
# verify_docs.py
def verify_parameter(param_name, doc_value, file_pattern):
    """Check if parameter exists in code with documented default"""
    # Search for parameter in source
    # Extract default value
    # Compare with documentation
    # Return pass/fail
```

### Phase 6: Graphics and Visualization (Day 8)

#### 6.1 Fix All Flowcharts
- [ ] Install mermaid-cli locally
- [ ] Convert all flowcharts to SVG
- [ ] Create clickable image maps
- [ ] Add alt text for accessibility
- [ ] Test in multiple browsers

#### 6.2 Create New Diagrams
- [ ] Module dependency graph
- [ ] Data structure hierarchy
- [ ] Task execution timeline
- [ ] AMR tree visualization
- [ ] Parallel decomposition diagram

### Phase 7: Final Review and Testing (Day 9-10)

#### 7.1 Comprehensive Review Checklist
- [ ] Every source file has documentation
- [ ] Every parameter verified against code
- [ ] All examples tested and working
- [ ] All flowcharts rendering correctly
- [ ] No HDF5 references remain
- [ ] Binary output properly documented
- [ ] All 20+ modules documented
- [ ] Cross-references validated

#### 7.2 Team Testing Protocol
```markdown
## Documentation Test for Team Members

1. Pick a random source file
2. Find its documentation
3. Verify accuracy against code
4. Try to build/run an example
5. Report discrepancies

Target: 100% accuracy
```

## Success Metrics

### Must Have (Critical)
- ✅ All flowcharts render correctly
- ✅ No factual errors (HDF5 removed)
- ✅ All 20+ modules documented
- ✅ Correct output formats listed
- ✅ Accurate parameter defaults

### Should Have (Important)
- ✅ Every source file referenced
- ✅ Working examples for each module
- ✅ Dependency graphs accurate
- ✅ Build instructions verified

### Nice to Have
- ✅ Interactive diagrams
- ✅ Video tutorials
- ✅ Jupyter notebooks

## Input Parameter Verification Script

```python
#!/usr/bin/env python3
# verify_input_params.py

import re
import os
from pathlib import Path

def scan_for_parameters(src_dir):
    """Scan source code for all parameter reads"""
    params = {}
    
    # Patterns to match parameter reading
    patterns = [
        r'GetString\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)',
        r'GetInteger\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)',
        r'GetReal\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)',
        r'GetBoolean\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)',
        r'GetOrAdd\w+\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*([^)]+)\)',
    ]
    
    for cpp_file in Path(src_dir).rglob('*.cpp'):
        with open(cpp_file, 'r') as f:
            content = f.read()
            for pattern in patterns:
                matches = re.findall(pattern, content)
                for match in matches:
                    block = match[0]
                    param = match[1]
                    if block not in params:
                        params[block] = {}
                    params[block][param] = {
                        'file': str(cpp_file),
                        'default': match[2] if len(match) > 2 else None
                    }
    
    return params

def generate_param_docs(params):
    """Generate markdown documentation for all parameters"""
    for block, block_params in sorted(params.items()):
        print(f"\n## Input Block: <{block}>")
        print("\n| Parameter | Default | Source File |")
        print("|-----------|---------|-------------|")
        for param, info in sorted(block_params.items()):
            default = info['default'] or 'required'
            file = info['file'].split('/')[-1]
            print(f"| {param} | {default} | {file} |")
```

## Execution Timeline

| Day | Task | Verification Method |
|-----|------|-------------------|
| 1 | Fix flowcharts, create master diagrams | Visual inspection |
| 2 | Document ALL input parameters | Script verification |
| 3 | Systematic content verification | Checklist completion |
| 4-5 | Complete all module documentation | Source code review |
| 6 | File-by-file reference | Automated script |
| 7 | Testing all examples | Run each example |
| 8 | Fix graphics and visualizations | Browser testing |
| 9-10 | Final review and team testing | Team feedback |

## Definition: What is a "Module" in AthenaK?

A **module** in AthenaK is a self-contained component that:
1. Has its own directory in `src/`
2. Implements a specific physics solver or numerical method
3. Can be enabled/disabled at compile time or runtime
4. Has its own task list entries
5. Manages its own data structures

**Not modules** but important components:
- Riemann solvers (subfunctions of physics modules)
- Reconstruction methods (numerical tools used by modules)
- Source terms (additions to physics modules)

## Quality Assurance Process

### For Each Documentation Page:
1. **Write**: Based on source code, not assumptions
2. **Verify**: Check against actual implementation
3. **Test**: Run examples to confirm behavior
4. **Review**: Have another person verify
5. **Update**: Fix any discrepancies found

### Documentation Rules:
- **NEVER** document features that don't exist
- **ALWAYS** check source code first
- **VERIFY** every parameter and default value
- **TEST** every example before including
- **CITE** source file and line numbers

## Deliverables

1. **Fixed flowcharts** - All rendering correctly as SVG/PNG
2. **20+ module docs** - Complete, accurate, tested
3. **File reference** - Every source file documented
4. **Corrected outputs** - No HDF5, proper binary/VTK docs
5. **Verification log** - Proof of testing
6. **Team guide** - How to maintain docs

## Master Verification Checklist

### Critical Fixes (Must Complete)
- [ ] All Mermaid flowcharts render as visual diagrams, not text
- [ ] Remove ALL references to HDF5 (does not exist)
- [ ] Document binary output format as primary option
- [ ] Add `single_file_per_rank` parameter documentation
- [ ] Document all 14 actual output formats (vtk, bin, rst, tab, cbin, pvtk, trk, df, dxh, ppd, prst, hst, log, pdf)

### Module Documentation (20 Total Required)
- [ ] Core: mesh, driver, tasklist, coordinates (4)
- [ ] Physics: hydro, mhd, radiation, z4c, dyn_grmhd, ion-neutral, particles (7)
- [ ] Numerical: reconstruction, riemann_solvers, eos, diffusion (4)
- [ ] Support: outputs, boundaries, srcterms, shearing_box, pgen (5)

### Input Parameter Documentation
- [ ] Every `<block>` documented with source file reference
- [ ] Every parameter with type, default, and line number
- [ ] Cross-references to modules that use each block
- [ ] Example .athinput file for each major use case

### File-by-File Reference
- [ ] Every .cpp file has description
- [ ] Every .hpp file has description
- [ ] Dependencies mapped
- [ ] Public functions listed

### Flowchart Requirements
- [ ] Main execution flow (main.cpp → driver → output)
- [ ] Task dependency graph
- [ ] AMR refinement process
- [ ] Module interaction diagram
- [ ] Data structure hierarchy
- [ ] All flowcharts clickable with navigation

### Testing Requirements
- [ ] Every documented example builds successfully
- [ ] Every documented parameter works as described
- [ ] Every output format tested and verified
- [ ] Navigation links all functional
- [ ] Search functionality works

## Final Validation Protocol

```bash
#!/bin/bash
# validate_docs.sh

echo "=== AthenaK Documentation Validation ==="

# 1. Check flowcharts render
echo "Checking flowcharts..."
for file in docs/build/html/flowcharts/*.html; do
    if grep -q "mermaid.initialize" "$file"; then
        echo "✓ Flowchart in $file"
    else
        echo "✗ FAILED: $file missing flowchart"
    fi
done

# 2. Check for HDF5 references
echo "Checking for incorrect HDF5 references..."
if grep -r "hdf5\|HDF5" docs/source/ --include="*.md" --include="*.rst"; then
    echo "✗ FAILED: HDF5 references found (must remove)"
else
    echo "✓ No HDF5 references"
fi

# 3. Verify output formats
echo "Checking output format documentation..."
for format in vtk bin rst tab cbin pvtk trk df dxh ppd prst hst log pdf; do
    if grep -q "$format" docs/source/modules/outputs.md; then
        echo "✓ Format $format documented"
    else
        echo "✗ MISSING: Format $format"
    fi
done

# 4. Count modules
MODULE_COUNT=$(find docs/source/modules -name "*.md" | wc -l)
echo "Module documentation: $MODULE_COUNT/20"

# 5. Test examples
echo "Testing documented examples..."
# Build and run each example

echo "=== Validation Complete ==="
```

## Next Steps

1. **Immediate**: Fix Mermaid rendering issue using the detailed fix process
2. **Hour 2**: Remove all HDF5 references, document actual output formats
3. **Hour 3**: Run parameter scanning script to document all input blocks
4. **Day 2**: Complete all 20 module documentations
5. **Day 3**: Run full validation protocol
6. **Day 4**: Team review with validation checklist

This plan ensures comprehensive, accurate documentation that reflects the actual codebase, not theoretical capabilities. Every piece will be verified against source code.