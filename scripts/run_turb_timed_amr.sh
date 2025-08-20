#!/bin/bash
#
# Script to run all 6 stages of the time-dependent AMR turbulence simulation
# Each stage refines at specific times and then stops for restart with adjusted parameters
#
# Usage: ./run_turb_timed_amr.sh [build_dir]
#
# Default build directory is ../build

# Set build directory (default to ../build if not specified)
BUILD_DIR=${1:-"../build"}
ATHENA_EXEC="${BUILD_DIR}/src/athena"
INPUT_DIR="../inputs/dbf"

# Check if athena executable exists
if [ ! -f "${ATHENA_EXEC}" ]; then
    echo "Error: athena executable not found at ${ATHENA_EXEC}"
    echo "Please build AthenaK with -DPROBLEM=turb_timed_amr first"
    exit 1
fi

# Function to run a stage
run_stage() {
    local stage=$1
    local cores=$2
    local input_file="${INPUT_DIR}/turb_timed_amr_stage${stage}.athinput"
    
    echo "=================================================="
    echo "Running Stage ${stage} with ${cores} cores"
    echo "Input file: ${input_file}"
    echo "=================================================="
    
    if [ ! -f "${input_file}" ]; then
        echo "Error: Input file ${input_file} not found"
        exit 1
    fi
    
    # Check if this is a restart (stage > 1)
    if [ $stage -gt 1 ]; then
        prev_stage=$((stage - 1))
        restart_file="TurbTimedAMR_S${prev_stage}.00000.rst"
        if [ ! -f "${restart_file}" ]; then
            echo "Warning: Restart file ${restart_file} not found"
            echo "Make sure stage ${prev_stage} completed successfully"
        fi
    fi
    
    # Run with MPI if cores > 1
    if [ $cores -gt 1 ]; then
        echo "Running with MPI: mpirun -np ${cores} ${ATHENA_EXEC} -i ${input_file}"
        mpirun -np ${cores} ${ATHENA_EXEC} -i ${input_file}
    else
        echo "Running serial: ${ATHENA_EXEC} -i ${input_file}"
        ${ATHENA_EXEC} -i ${input_file}
    fi
    
    # Check if run completed successfully
    if [ $? -eq 0 ]; then
        echo "Stage ${stage} completed successfully"
    else
        echo "Error: Stage ${stage} failed"
        exit 1
    fi
    
    echo ""
}

# Main execution
echo "Starting Time-Dependent AMR Turbulence Simulation"
echo "=================================================="
echo ""

# Stage 1: Base resolution (64^3 -> 128^3)
echo "Stage 1: Initial run to t=2.0 eddy times"
echo "Resolution: 64^3 base, will refine to 128^3"
echo "Suggested cores: 8"
run_stage 1 8

# Stage 2: First refinement (128^3 -> 256^3)
echo "Stage 2: Continue to t=3.0 eddy times"
echo "Resolution: 128^3 base, will refine to 256^3"
echo "Suggested cores: 16"
echo "Note: Viscosity reduced by factor of 2"
run_stage 2 16

# Stage 3: Second refinement (256^3 -> 512^3)
echo "Stage 3: Continue to t=3.5 eddy times"
echo "Resolution: 256^3 base, will refine to 512^3"
echo "Suggested cores: 32"
echo "Note: Viscosity reduced by factor of 4"
run_stage 3 32

# Stage 4: Third refinement (512^3 -> 1024^3)
echo "Stage 4: Continue to t=3.75 eddy times"
echo "Resolution: 512^3 base, will refine to 1024^3"
echo "Suggested cores: 64"
echo "Note: Viscosity reduced by factor of 8"
run_stage 4 64

# Stage 5: Fourth refinement (1024^3 -> 2048^3)
echo "Stage 5: Continue to t=3.875 eddy times"
echo "Resolution: 1024^3 base, will refine to 2048^3"
echo "Suggested cores: 128"
echo "Note: Viscosity reduced by factor of 16"
run_stage 5 128

# Stage 6: Final run at highest resolution
echo "Stage 6: Final run to t=5.875 eddy times"
echo "Resolution: 2048^3 (no further refinement)"
echo "Suggested cores: 128-256"
echo "Note: Running 2 eddy times at highest resolution"
run_stage 6 256

echo "=================================================="
echo "All stages completed successfully!"
echo "Final output files: TurbTimedAMR_S6.*.vtk"
echo "=================================================="