#!/bin/bash
# Example usage of the SFB turbulence analysis script

# Set paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ANALYSIS_SCRIPT="${SCRIPT_DIR}/analyze_sfb_output.py"
OUTPUT_DIR="${SCRIPT_DIR}/sfb_analysis_output"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "SFB Turbulence Analysis Example"
echo "==============================="
echo ""
echo "This script demonstrates how to analyze SFB turbulence simulation output."
echo ""

# Check if we have any output files
if [ ! -f "*.athdf" ] && [ ! -f "*.hdf5" ]; then
    echo "No AthenaK output files found in current directory."
    echo ""
    echo "To generate output files, run:"
    echo "  cd /path/to/athena/build"
    echo "  ./src/athena -i inputs/hydro/sfb_turb.athinput"
    echo ""
    echo "Then run this script in the directory containing the output files."
    exit 1
fi

# Single file analysis
echo "1. Analyzing single output file..."
if [ -f "sfb_turb.out1.00000.athdf" ]; then
    python3 "${ANALYSIS_SCRIPT}" sfb_turb.out1.00000.athdf \
        --r0_turb 0.4 \
        --output-dir "${OUTPUT_DIR}" \
        --power-spectrum
fi

# Time series analysis with history file
echo ""
echo "2. Analyzing time series..."
if [ -f "sfb_turb.hst" ]; then
    python3 "${ANALYSIS_SCRIPT}" "sfb_turb.out1.*.athdf" \
        --r0_turb 0.4 \
        --history sfb_turb.hst \
        --output-dir "${OUTPUT_DIR}" \
        --all-frames
fi

# Power spectrum analysis on final frame
echo ""
echo "3. Computing power spectrum from final frame..."
LAST_FILE=$(ls -1 sfb_turb.out1.*.athdf 2>/dev/null | tail -1)
if [ -n "${LAST_FILE}" ]; then
    python3 "${ANALYSIS_SCRIPT}" "${LAST_FILE}" \
        --r0_turb 0.4 \
        --power-spectrum \
        --expo 1.667 \
        --output-dir "${OUTPUT_DIR}"
fi

echo ""
echo "Analysis complete! Results saved to: ${OUTPUT_DIR}"
echo ""
echo "Generated plots:"
echo "  - *_slice.png: 2D slices showing density, velocity, vorticity, etc."
echo "  - *_radial.png: Spherically averaged radial profiles"
echo "  - power_spectrum.png: Velocity power spectrum"
echo "  - time_series.png: Time evolution of key quantities"

# Additional analysis suggestions
echo ""
echo "Additional analysis options:"
echo "  --r0_turb VALUE     Set turbulence region radius (default: 0.4)"
echo "  --expo VALUE        Expected power spectrum slope (default: 5/3)"
echo "  --all-frames        Analyze all frames (for animations)"
echo ""
echo "Example custom analysis:"
echo "  python3 ${ANALYSIS_SCRIPT} output.athdf --r0_turb 0.3 --expo 2.0"