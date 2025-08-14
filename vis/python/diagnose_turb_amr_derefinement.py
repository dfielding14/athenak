#!/usr/bin/env python3
"""
Diagnostic script to analyze turbulence driving behavior during AMR derefinement.
This script tracks:
1. MeshBlock count changes over time
2. Force field continuity at refinement boundaries
3. Array indexing issues during derefinement
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import athena_read as ar
import glob
import sys

def read_binary_as_athdf(filename):
    """Read AthenaK binary output file"""
    return ar.athdf(filename)

def analyze_meshblock_evolution(output_dir):
    """Track MeshBlock count and structure over time"""
    
    # Get all output files
    mb_files = sorted(glob.glob(f"{output_dir}/TurbAMRWave.out3.*.bin"))
    
    if len(mb_files) == 0:
        print("No MeshBlock output files found!")
        return None, None, None
    
    times = []
    mb_counts = []
    mb_levels = []
    
    for file in mb_files:
        try:
            data = read_binary_as_athdf(file)
            time = float(data['Time'])
            
            # Get MeshBlock GIDs
            mb_gid = data['mb_gid']
            
            # Count total MeshBlocks
            nmb = len(np.unique(mb_gid))
            
            # Estimate refinement levels from MeshBlock sizes
            # Smaller blocks indicate higher refinement
            x1f = data['x1f']
            dx1_per_mb = []
            
            for mb in range(data['MeshBlockSize'][0]):
                i_start = mb * data['MeshBlockSize'][1]
                i_end = (mb + 1) * data['MeshBlockSize'][1]
                if i_end <= len(x1f):
                    dx1 = x1f[i_end-1] - x1f[i_start]
                    dx1_per_mb.append(dx1)
            
            times.append(time)
            mb_counts.append(nmb)
            mb_levels.append(dx1_per_mb)
            
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    
    return np.array(times), np.array(mb_counts), mb_levels

def analyze_force_continuity(output_dir, time_index):
    """Check force field continuity at a specific time"""
    
    force_files = sorted(glob.glob(f"{output_dir}/TurbAMRWave.out2.*.bin"))
    
    if time_index >= len(force_files):
        print(f"Time index {time_index} out of range")
        return None
    
    try:
        data = read_binary_as_athdf(force_files[time_index])
        
        # Get force components
        force1 = data['force1']
        force2 = data['force2']
        force3 = data['force3']
        
        # Calculate force magnitude
        force_mag = np.sqrt(force1**2 + force2**2 + force3**2)
        
        # Get coordinates
        x1v = data['x1v']
        x2v = data['x2v']
        x3v = data['x3v']
        
        return {
            'time': float(data['Time']),
            'x1v': x1v,
            'x2v': x2v,
            'x3v': x3v,
            'force_mag': force_mag,
            'force1': force1,
            'force2': force2,
            'force3': force3
        }
        
    except Exception as e:
        print(f"Error analyzing force continuity: {e}")
        return None

def plot_diagnostics(times, mb_counts, force_data_list):
    """Create diagnostic plots"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: MeshBlock count over time
    ax = axes[0, 0]
    if times is not None and mb_counts is not None:
        ax.plot(times, mb_counts, 'b.-', linewidth=2, markersize=8)
        ax.set_xlabel('Time')
        ax.set_ylabel('Number of MeshBlocks')
        ax.set_title('MeshBlock Count Evolution')
        ax.grid(True, alpha=0.3)
        
        # Highlight derefinement events (where MB count decreases)
        for i in range(1, len(mb_counts)):
            if mb_counts[i] < mb_counts[i-1]:
                ax.axvline(times[i], color='red', alpha=0.5, linestyle='--',
                          label='Derefinement' if i == 1 else '')
        if 'Derefinement' in [l.get_label() for l in ax.lines]:
            ax.legend()
    
    # Plot 2: Force magnitude at different times
    ax = axes[0, 1]
    colors = plt.cm.viridis(np.linspace(0, 1, len(force_data_list)))
    
    for i, force_data in enumerate(force_data_list):
        if force_data is not None:
            # Take a slice through y=0, z=0
            ny = len(force_data['x2v'])
            nz = len(force_data['x3v'])
            j_mid = ny // 2
            k_mid = nz // 2
            
            force_slice = force_data['force_mag'][k_mid, j_mid, :]
            
            ax.plot(force_data['x1v'], force_slice, 
                   color=colors[i], label=f"t={force_data['time']:.2f}")
    
    ax.set_xlabel('x1')
    ax.set_ylabel('Force Magnitude')
    ax.set_title('Force Profile Evolution (y=0, z=0 slice)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: 2D force magnitude at a specific time
    if len(force_data_list) > 0 and force_data_list[-1] is not None:
        ax = axes[1, 0]
        force_data = force_data_list[-1]
        
        # Take z=0 slice
        nz = len(force_data['x3v'])
        k_mid = nz // 2
        
        X1, X2 = np.meshgrid(force_data['x1v'], force_data['x2v'])
        force_2d = force_data['force_mag'][k_mid, :, :]
        
        im = ax.pcolormesh(X1, X2, force_2d, shading='nearest', cmap='RdBu_r')
        ax.set_xlabel('x1')
        ax.set_ylabel('x2')
        ax.set_title(f'Force Magnitude at t={force_data["time"]:.2f} (z=0 slice)')
        ax.set_aspect('equal')
        plt.colorbar(im, ax=ax, label='|F|')
    
    # Plot 4: Force discontinuity metric
    ax = axes[1, 1]
    
    discontinuity_times = []
    discontinuity_values = []
    
    for force_data in force_data_list:
        if force_data is not None:
            # Calculate gradient magnitude as discontinuity metric
            force_mag = force_data['force_mag']
            
            # Take central differences in x1 direction
            grad_x1 = np.zeros_like(force_mag)
            grad_x1[:, :, 1:-1] = (force_mag[:, :, 2:] - force_mag[:, :, :-2]) / 2.0
            
            # Maximum gradient as discontinuity metric
            max_grad = np.max(np.abs(grad_x1))
            
            discontinuity_times.append(force_data['time'])
            discontinuity_values.append(max_grad)
    
    if len(discontinuity_times) > 0:
        ax.plot(discontinuity_times, discontinuity_values, 'r.-', linewidth=2)
        ax.set_xlabel('Time')
        ax.set_ylabel('Max Force Gradient')
        ax.set_title('Force Discontinuity Metric')
        ax.grid(True, alpha=0.3)
        
        # Mark derefinement times
        if times is not None and mb_counts is not None:
            for i in range(1, len(mb_counts)):
                if mb_counts[i] < mb_counts[i-1]:
                    ax.axvline(times[i], color='red', alpha=0.5, linestyle='--')
    
    plt.tight_layout()
    return fig

def main():
    if len(sys.argv) < 2:
        print("Usage: python diagnose_turb_amr_derefinement.py <output_directory>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    
    print("Analyzing turbulence AMR derefinement diagnostics...")
    
    # Analyze MeshBlock evolution
    times, mb_counts, mb_levels = analyze_meshblock_evolution(output_dir)
    
    if times is None:
        print("Failed to analyze MeshBlock evolution")
        return
    
    print(f"Time range: {times[0]:.2f} to {times[-1]:.2f}")
    print(f"MeshBlock count range: {np.min(mb_counts)} to {np.max(mb_counts)}")
    
    # Find derefinement events
    derefine_times = []
    for i in range(1, len(mb_counts)):
        if mb_counts[i] < mb_counts[i-1]:
            derefine_times.append(times[i])
            print(f"Derefinement at t={times[i]:.2f}: {mb_counts[i-1]} -> {mb_counts[i]} MeshBlocks")
    
    # Analyze force field at several time points, focusing on derefinement times
    force_data_list = []
    
    # Sample times: beginning, middle, derefinement events, and end
    sample_indices = [0]  # Start
    
    # Add derefinement time indices
    for dt in derefine_times:
        idx = np.argmin(np.abs(times - dt))
        if idx not in sample_indices:
            sample_indices.append(idx)
    
    # Add a few more samples
    sample_indices.extend([len(times)//3, 2*len(times)//3, len(times)-1])
    sample_indices = sorted(list(set(sample_indices)))[:6]  # Limit to 6 samples
    
    print("\nAnalyzing force field at times:")
    for idx in sample_indices:
        if idx < len(times):
            print(f"  t={times[idx]:.2f} (index {idx})")
            force_data = analyze_force_continuity(output_dir, idx)
            force_data_list.append(force_data)
    
    # Create diagnostic plots
    fig = plot_diagnostics(times, mb_counts, force_data_list)
    
    output_file = f"{output_dir}/turb_amr_derefinement_diagnostics.png"
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nDiagnostic plots saved to: {output_file}")
    
    # Additional analysis: Check for array size mismatches
    print("\nPotential issues to investigate:")
    print("1. Turbulence driver arrays may not resize during mesh refinement")
    print("2. MeshBlock indices in force arrays may become invalid after derefinement")
    print("3. Phase information may be lost when MeshBlocks are destroyed/created")
    
    if len(derefine_times) > 0:
        print("\nRecommendation: Add debug output in TurbulenceDriver to track:")
        print("  - Array sizes vs actual MeshBlock count")
        print("  - MeshBlock index mapping before/after refinement")
        print("  - Force field values at refinement boundaries")

if __name__ == "__main__":
    main()