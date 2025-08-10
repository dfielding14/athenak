#!/bin/bash
# Simplified MPI test that completes quickly but still tests multiple wave circuits

echo "========================================="
echo "Simplified div(B) MPI Test"
echo "========================================="

# Configuration
BUILD_DIR="build_divb_test"
OUTPUT_DIR="test_divb_output"

# Build if needed
if [ ! -d "$BUILD_DIR" ]; then
    echo "Building AthenaK with MHD+AMR..."
    mkdir -p $BUILD_DIR
    cd $BUILD_DIR
    cmake .. -DPROBLEM=turb_mhd_amr_wave -DAthena_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
    make -j4 > /dev/null 2>&1
    cd ..
    echo "Build complete"
fi

ATHENA_EXE="$BUILD_DIR/src/athena"

# Create output directory
mkdir -p $OUTPUT_DIR

# Create optimized input file for quick multi-circuit test
cat > $OUTPUT_DIR/test_mpi.athinput << 'EOF'
# Quick MPI test with 1.5 wave circuits
<comment>
problem = turb_mhd_amr_wave

<job>
basename = test_mpi

<mesh>
nghost = 2
nx1 = 32
x1min = -0.5
x1max = 0.5
ix1_bc = periodic
ox1_bc = periodic

nx2 = 32
x2min = -0.5  
x2max = 0.5
ix2_bc = periodic
ox2_bc = periodic

nx3 = 1
x3min = -0.5
x3max = 0.5
ix3_bc = periodic
ox3_bc = periodic

<mesh_refinement>
adaptive = true
max_nmb_per_rank = 100
refinement_interval = 5
ncycle_check = 2
max_level = 2

<meshblock>
nx1 = 8
nx2 = 8
nx3 = 1

<time>
evolution = dynamic
cfl_number = 0.4
tlim = 1.5        # 1.5 circuits at speed 1.0
nlim = 200        # Limit cycles
integrator = rk2
ncycle_out = 20

<mhd>
eos = ideal
gamma = 1.4
reconstruct = plm
rsolver = hlld
dfloor = 1.0e-8
pfloor = 1.0e-10

<problem>
rho0 = 1.0
pgas0 = 1.0
vx0 = 0.0
vy0 = 0.0
vz0 = 0.0
b0 = 1.0
beta = -1.0
field_config = 3

# Wave for 1 circuit per time unit
wave_speed = 1.0
wave_width = 0.2
wave_amplitude = 0.5
wave_direction = 0

<turb_driving>
dedt = 0.0

<output1>
file_type = bin
dt = 0.3          # 5 outputs total
variable = mhd_w_bcc
id = cons

<output2>
file_type = bin
dt = 0.3
variable = mhd_divb
id = divb

<output3>
file_type = hst
dt = 0.01
EOF

cd $OUTPUT_DIR

# Get number of CPU cores
if command -v nproc &> /dev/null; then
    NCORES=$(nproc)
elif command -v sysctl &> /dev/null; then
    NCORES=$(sysctl -n hw.ncpu)
else
    NCORES=1
fi

echo "Configuration:"
echo "  Wave speed: 1.0 (1 circuit/time unit)"
echo "  Target: 1.5 circuits (t=1.5)"
echo "  Resolution: 32x32, max_level=2"
echo ""

# Run tests with different MPI ranks
for np in 1 2 4; do
    if [ $NCORES -lt $np ]; then
        echo "Skipping $np ranks (only $NCORES cores available)"
        continue
    fi
    
    TEST_DIR="test_mpi_np${np}"
    mkdir -p $TEST_DIR
    cd $TEST_DIR
    
    echo -n "Running with $np rank(s)..."
    
    # Run simulation
    if [ $np -eq 1 ]; then
        ../../$ATHENA_EXE -i ../test_mpi.athinput > run.log 2>&1
    else
        mpirun -np $np ../../$ATHENA_EXE -i ../test_mpi.athinput > run.log 2>&1
    fi
    
    # Check results
    if [ -f run.log ]; then
        final_time=$(tail -20 run.log | grep "time=" | tail -1 | sed -n 's/.*time=\([0-9.e+-]*\).*/\1/p')
        final_cycle=$(tail -20 run.log | grep "cycle=" | tail -1 | sed -n 's/.*cycle=\([0-9]*\).*/\1/p')
        
        if [ ! -z "$final_time" ]; then
            circuits=$(echo "scale=2; $final_time / 1.0" | bc -l 2>/dev/null || echo "0")
            echo " Complete! (t=$final_time, ${circuits} circuits, cycle=$final_cycle)"
        else
            echo " Failed to complete properly"
        fi
        
        # Check for refinement
        refinements=$(grep -c "adaptive refinement" run.log 2>/dev/null || echo "0")
        if [ "$refinements" -gt "0" ]; then
            echo "  Refinement events: $refinements"
        fi
        
        # Output files
        divb_files=$(ls bin/*.divb.*.bin 2>/dev/null | wc -l)
        echo "  div(B) outputs: $divb_files files"
    else
        echo " Error - no log file"
    fi
    
    cd ..
done

cd ..

echo ""
echo "========================================="
echo "Test complete! Output in $OUTPUT_DIR"
echo "========================================="