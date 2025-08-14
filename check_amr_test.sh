#!/bin/bash
cd /Users/dbf75/Work/Research/athenak-DF/build_wave/src

# Run the test and capture output for analysis
echo "Running AMR test..."
timeout 30s ./athena -i ../../inputs/hydro/turb_amr_wave_test.athinput 2>&1 | tee amr_test_output.log

# Check for refinement activity
echo -e "\n=== Refinement Summary ==="
grep -E "(MeshBlocks created|deleted by AMR)" amr_test_output.log | tail -5
echo -e "\n=== Refinement Conditions ==="
grep -E "(MB [0-9]:.*level=)" amr_test_output.log | tail -10
echo -e "\n=== Resize Activity ==="
grep -E "(CheckResize|ResizeArrays)" amr_test_output.log | tail -5

# Clean up
rm -f amr_test_output.log