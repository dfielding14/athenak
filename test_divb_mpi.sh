#!/bin/bash
#==============================================================================
# test_divb_mpi.sh
# Test script for verifying div(B) preservation with multiple MPI ranks
# This script runs the MHD+AMR tests with different numbers of MPI processes
# and compares the results to check for cross-rank boundary issues
#==============================================================================

# Configuration
ATHENA_EXE="./athena"
INPUT_DIR="inputs/mhd"
OUTPUT_DIR="test_divb_output"
BUILD_DIR="build_test"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "========================================="
echo "div(B) AMR Multi-Rank Test Suite"
echo "========================================="

# Check if executable exists
if [ ! -f "$BUILD_DIR/src/athena" ]; then
    echo -e "${RED}Error: Athena executable not found at $BUILD_DIR/src/athena${NC}"
    echo "Building Athena with MHD+AMR test problem..."
    
    # Build with MPI support
    mkdir -p $BUILD_DIR
    cd $BUILD_DIR
    cmake .. -DPROBLEM=turb_mhd_amr_wave -DAthena_ENABLE_MPI=ON
    make -j4
    cd ..
    
    if [ ! -f "$BUILD_DIR/src/athena" ]; then
        echo -e "${RED}Build failed!${NC}"
        exit 1
    fi
fi

ATHENA_EXE="$BUILD_DIR/src/athena"

# Create output directory
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# Function to run test with given number of ranks
run_test() {
    local nprocs=$1
    local test_name=$2
    local input_file=$3
    
    echo ""
    echo -e "${GREEN}Running $test_name with $nprocs MPI rank(s)...${NC}"
    
    # Create subdirectory for this test
    local test_dir="${test_name}_np${nprocs}"
    mkdir -p $test_dir
    cd $test_dir
    
    # Run the simulation
    if [ $nprocs -eq 1 ]; then
        ../../$ATHENA_EXE -i ../../$input_file 2>&1 | tee run.log
    else
        mpirun -np $nprocs ../../$ATHENA_EXE -i ../../$input_file 2>&1 | tee run.log
    fi
    
    # Check if run completed successfully
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Test completed successfully${NC}"
        
        # Extract div(B) diagnostics from log
        echo ""
        echo "div(B) Fix Diagnostics:"
        grep "\[FixDivB\]" run.log | tail -20
        
        # Count fixed boundaries
        local total_fixed=$(grep "\[FixDivB\]" run.log | \
                           awk '{sum += $4} END {print sum}')
        local total_cross=$(grep "\[FixDivB\]" run.log | \
                           awk '{sum += $7} END {print sum}')
        
        echo ""
        echo "Summary: Fixed $total_fixed same-rank boundaries"
        if [ $nprocs -gt 1 ]; then
            echo "         Detected $total_cross potential cross-rank boundaries"
        fi
    else
        echo -e "${RED}✗ Test failed!${NC}"
    fi
    
    cd ..
}

# Test 1: 2D test with different numbers of ranks
echo ""
echo -e "${YELLOW}Test Set 1: 2D MHD+AMR with traveling wave${NC}"
echo "==========================================="

for np in 1 2 4 8; do
    # Check if we have enough cores
    if [ $(nproc) -ge $np ]; then
        run_test $np "2d_wave" "$INPUT_DIR/turb_mhd_amr_wave_2d.athinput"
    else
        echo -e "${YELLOW}Skipping $np ranks (only $(nproc) cores available)${NC}"
    fi
done

# Test 2: 1D test (quick test)
echo ""
echo -e "${YELLOW}Test Set 2: 1D MHD+AMR (quick test)${NC}"
echo "====================================="

for np in 1 2 4; do
    if [ $(nproc) -ge $np ]; then
        run_test $np "1d_wave" "$INPUT_DIR/turb_mhd_amr_wave_1d.athinput"
    fi
done

# Test 3: 3D test (if enough memory/cores)
echo ""
echo -e "${YELLOW}Test Set 3: 3D MHD+AMR (optional, resource intensive)${NC}"
echo "====================================================="
echo "Run manually with: mpirun -np N $ATHENA_EXE -i $INPUT_DIR/turb_mhd_amr_wave.athinput"

# Summary
echo ""
echo "========================================="
echo "Test Suite Complete"
echo "========================================="
echo "Output files are in: $OUTPUT_DIR"
echo ""
echo "To analyze div(B) in output files, run:"
echo "  python3 analyze_divb.py $OUTPUT_DIR"
echo ""
echo "Key things to look for:"
echo "1. Compare total fixed boundaries between 1 rank and N ranks"
echo "2. Check if potential cross-rank boundaries increase with more ranks"
echo "3. Examine div(B) values in output HDF5 files near rank boundaries"