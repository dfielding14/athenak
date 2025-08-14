#!/bin/bash
# Comprehensive validation of AthenaK documentation

echo "========================================="
echo "   AthenaK Documentation Full Validation"
echo "========================================="

FAILED=0
PASSED=0

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

check() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $2"
        ((PASSED++))
    else
        echo -e "${RED}✗${NC} $2"
        ((FAILED++))
    fi
}

echo
echo "1. Checking for HDF5 references (should be NONE)..."
HDF5_COUNT=$(grep -r "hdf5\|HDF5" source/ --include="*.md" --include="*.rst" 2>/dev/null | wc -l)
if [ $HDF5_COUNT -eq 0 ]; then
    check 0 "No HDF5 references found (correct)"
else
    check 1 "Found $HDF5_COUNT HDF5 references (MUST REMOVE)"
    grep -r "hdf5\|HDF5" source/ --include="*.md" --include="*.rst"
fi

echo
echo "2. Verifying output formats documentation..."
for format in vtk bin rst tab cbin pvtk trk df dxh ppd prst hst log pdf; do
    if grep -q "$format" source/modules/outputs.md 2>/dev/null; then
        check 0 "Format '$format' documented"
    else
        check 1 "Format '$format' MISSING"
    fi
done

echo
echo "3. Checking single_file_per_rank documentation..."
if grep -q "single_file_per_rank" source/modules/outputs.md 2>/dev/null; then
    check 0 "single_file_per_rank documented"
else
    check 1 "single_file_per_rank NOT documented"
fi

echo
echo "4. Verifying module documentation..."
EXPECTED_MODULES=(mesh driver tasklist coordinates hydro mhd radiation z4c dyn_grmhd 
                  ion_neutral particles reconstruction riemann_solvers eos diffusion 
                  outputs boundaries srcterms shearing_box pgen)

for module in "${EXPECTED_MODULES[@]}"; do
    if [ -f "source/modules/${module}.md" ]; then
        check 0 "Module '$module' documented"
    else
        check 1 "Module '$module' MISSING"
    fi
done

echo
echo "5. Checking flowchart rendering..."
if [ -f "build/html/flowcharts/runtime.html" ]; then
    # Check if mermaid scripts are included
    if grep -q "mermaid.initialize" build/html/flowcharts/runtime.html; then
        check 0 "Mermaid initialization present"
    else
        check 1 "Mermaid initialization MISSING"
    fi
    
    # Check for raw mermaid text (shouldn't be visible)
    if grep -q "flowchart TD" build/html/flowcharts/runtime.html | head -1; then
        echo -e "${YELLOW}⚠${NC} Flowchart may be rendering as text"
    fi
else
    check 1 "Runtime flowchart HTML not found"
fi

echo
echo "6. Verifying input parameters documentation..."
if [ -f "source/reference/input_parameters.md" ]; then
    BLOCK_COUNT=$(grep -c "^## Input Block:" source/reference/input_parameters.md)
    check 0 "Input parameters documented ($BLOCK_COUNT blocks)"
else
    check 1 "Input parameters documentation MISSING"
fi

echo
echo "7. Testing example configurations..."
# Check if documented examples are valid
echo "Testing example input blocks..."

# Test a basic hydro configuration
cat > test_hydro.athinput << EOF
<time>
integrator = rk2
cfl_number = 0.4
tlim = 1.0

<hydro>
eos = ideal
gamma = 1.4
reconstruct = plm
rsolver = hllc

<output1>
file_type = vtk
dt = 0.1
variable = prim
EOF

# Would need actual athena binary to test this
echo -e "${YELLOW}⚠${NC} Cannot test input files without athena binary"

echo
echo "8. Checking navigation links..."
BROKEN_LINKS=0
for file in build/html/*.html build/html/*/*.html; do
    if [ -f "$file" ]; then
        # Check for broken internal links
        grep -o 'href="[^"]*"' "$file" 2>/dev/null | while read -r link; do
            link_target=$(echo "$link" | sed 's/href="\([^"]*\)"/\1/')
            if [[ "$link_target" == *.html* ]] && [[ ! "$link_target" == http* ]]; then
                # Internal link - check if target exists
                target_file="build/html/$link_target"
                if [ ! -f "$target_file" ] && [[ ! "$link_target" == *#* ]]; then
                    ((BROKEN_LINKS++))
                fi
            fi
        done
    fi
done

if [ $BROKEN_LINKS -eq 0 ]; then
    check 0 "No broken internal links found"
else
    check 1 "Found $BROKEN_LINKS broken links"
fi

echo
echo "9. Verifying critical fixes..."
check_critical() {
    local desc="$1"
    local file="$2"
    local pattern="$3"
    local should_exist="$4"
    
    if [ "$should_exist" = "true" ]; then
        if grep -q "$pattern" "$file" 2>/dev/null; then
            check 0 "$desc"
        else
            check 1 "$desc NOT FOUND"
        fi
    else
        if grep -q "$pattern" "$file" 2>/dev/null; then
            check 1 "$desc STILL EXISTS (should be removed)"
        else
            check 0 "$desc"
        fi
    fi
}

check_critical "Binary output documented" "source/modules/outputs.md" "bin.*binary" true
check_critical "No HDF5 in outputs" "source/modules/outputs.md" "hdf5\|HDF5" false
check_critical "VTK format documented" "source/modules/outputs.md" "vtk.*ParaView" true

echo
echo "========================================="
echo "          VALIDATION SUMMARY"
echo "========================================="
echo -e "Passed: ${GREEN}$PASSED${NC}"
echo -e "Failed: ${RED}$FAILED${NC}"

if [ $FAILED -eq 0 ]; then
    echo -e "\n${GREEN}SUCCESS!${NC} Documentation passes all validations."
    exit 0
else
    echo -e "\n${RED}FAILURE!${NC} Documentation has $FAILED issues that need fixing."
    exit 1
fi