#!/bin/bash
#==============================================================================
# quick_test_divb.sh
# Quick test to check div(B) fix with multiple MPI ranks
# This runs a simple 2D test and immediately shows the diagnostic output
#==============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}======================================${NC}"
echo -e "${BLUE}Quick div(B) Cross-Rank Boundary Test${NC}"
echo -e "${BLUE}======================================${NC}"

# Step 1: Build
echo -e "\n${YELLOW}Step 1: Building AthenaK with MHD+AMR...${NC}"

BUILD_DIR="build_divb_test"
if [ ! -d "$BUILD_DIR" ]; then
    mkdir -p $BUILD_DIR
fi

cd $BUILD_DIR

# Configure with MPI
cmake .. -DPROBLEM=turb_mhd_amr_wave -DAthena_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Debug > /dev/null 2>&1

# Build
make -j4 > /dev/null 2>&1

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Build successful${NC}"
else
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi

cd ..

# Step 2: Run tests
echo -e "\n${YELLOW}Step 2: Running tests with different MPI ranks...${NC}"

ATHENA="$BUILD_DIR/src/athena"
INPUT="inputs/mhd/test_divb_ranks.athinput"
OUTPUT_BASE="test_quick"

# Function to run and analyze
run_and_analyze() {
    local np=$1
    local output_dir="${OUTPUT_BASE}_np${np}"
    
    echo -e "\n${GREEN}Running with $np MPI rank(s)...${NC}"
    
    # Clean previous output
    rm -rf $output_dir
    mkdir -p $output_dir
    cd $output_dir
    
    # Run simulation
    if [ $np -eq 1 ]; then
        ../$ATHENA -i ../$INPUT 2>&1 | tee run.log > /dev/null
    else
        mpirun -np $np ../$ATHENA -i ../$INPUT 2>&1 | tee run.log > /dev/null
    fi
    
    # Extract diagnostics
    echo -e "\n${BLUE}Diagnostic output:${NC}"
    grep "\[FixDivB\]" run.log | tail -10
    
    # Count boundaries
    if [ $np -gt 1 ]; then
        local total_fixed=$(grep "\[FixDivB\]" run.log 2>/dev/null | \
                           awk '{sum += $4} END {print sum}')
        local total_cross=$(grep "\[FixDivB\]" run.log 2>/dev/null | \
                           awk '{sum += $7} END {print sum}')
        
        echo -e "\n${YELLOW}Summary for $np ranks:${NC}"
        echo "  Total same-rank boundaries fixed: ${total_fixed:-0}"
        echo "  Potential cross-rank boundaries:  ${total_cross:-0}"
    else
        local total_fixed=$(grep "\[FixDivB\]" run.log 2>/dev/null | \
                           awk '{sum += $3} END {print sum}')
        echo -e "\n${YELLOW}Summary for 1 rank:${NC}"
        echo "  Total boundaries fixed: ${total_fixed:-0}"
    fi
    
    cd ..
}

# Get number of CPU cores (platform-aware)
if command -v nproc &> /dev/null; then
    NCORES=$(nproc)
elif command -v sysctl &> /dev/null; then
    NCORES=$(sysctl -n hw.ncpu)
else
    NCORES=1
fi

# Run with 1, 2, and 4 ranks
for np in 1 2 4; do
    if [ $NCORES -ge $np ]; then
        run_and_analyze $np
    else
        echo -e "${YELLOW}Skipping $np ranks (only $NCORES cores available)${NC}"
    fi
done

# Step 3: Analysis
echo -e "\n${YELLOW}Step 3: Analyzing div(B) in output files...${NC}"

# Quick Python analysis if available
if command -v python3 &> /dev/null; then
    python3 << 'EOF'
import os
import glob
import numpy as np

print("\nQuick div(B) analysis:")
print("-" * 40)

for np in [1, 2, 4]:
    test_dir = f"test_quick_np{np}"
    if not os.path.exists(test_dir):
        continue
    
    divb_files = glob.glob(f"{test_dir}/*divb*.athdf")
    if divb_files:
        # Would need h5py to read, just note existence
        print(f"np={np}: {len(divb_files)} div(B) output files created")
    
    # Check for any obvious errors in log
    log_file = f"{test_dir}/run.log"
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            log_content = f.read()
            if "ERROR" in log_content or "FATAL" in log_content:
                print(f"  ⚠️  Errors detected in run with {np} ranks!")

print("-" * 40)
EOF
else
    echo "Python3 not available for analysis"
fi

# Step 4: Summary
echo -e "\n${BLUE}======================================${NC}"
echo -e "${BLUE}Test Complete - Summary${NC}"
echo -e "${BLUE}======================================${NC}"

echo -e "\n${YELLOW}Key Observations:${NC}"
echo "1. Compare the number of fixed boundaries between 1 and N ranks"
echo "2. With N>1 ranks, check if 'potential cross-rank boundaries' > 0"
echo "3. If cross-rank boundaries exist, the fix is only partial"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "• To visualize div(B): python3 analyze_divb.py test_quick_np*"
echo "• To run longer test: use test_divb_mpi.sh"
echo "• To check specific rank boundaries: examine MeshBlock distribution"
echo ""
echo -e "${GREEN}Output directories:${NC}"
ls -d test_quick_np* 2>/dev/null || echo "No output directories found"