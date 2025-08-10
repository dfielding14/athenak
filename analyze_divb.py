#!/usr/bin/env python3
"""
analyze_divb.py
Analyze div(B) from AthenaK MHD+AMR simulation outputs
Compares div(B) errors between single-rank and multi-rank runs
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'vis', 'python'))
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
try:
    from bin_convert_new import read_binary, read_all_ranks_binary
except ImportError:
    print("Warning: Could not import bin_convert_new")
    print("Make sure vis/python/bin_convert_new.py is available")

def calculate_divb_from_fields(bx, by, bz, dx, dy, dz):
    """
    Calculate div(B) from face-centered magnetic fields
    
    Parameters:
    -----------
    bx, by, bz : arrays
        Face-centered magnetic field components
    dx, dy, dz : float
        Grid spacing in each direction
    
    Returns:
    --------
    divb : array
        Cell-centered div(B)
    """
    # bx has shape (nz, ny, nx+1) - face-centered in x
    # by has shape (nz, ny+1, nx) - face-centered in y  
    # bz has shape (nz+1, ny, nx) - face-centered in z
    
    # Calculate derivatives
    if bx.ndim == 3:  # 3D
        dbx_dx = (bx[:, :, 1:] - bx[:, :, :-1]) / dx
        dby_dy = (by[:, 1:, :] - by[:, :-1, :]) / dy
        dbz_dz = (bz[1:, :, :] - bz[:-1, :, :]) / dz
        divb = dbx_dx + dby_dy + dbz_dz
    elif bx.ndim == 2:  # 2D
        dbx_dx = (bx[:, 1:] - bx[:, :-1]) / dx
        dby_dy = (by[1:, :] - by[:-1, :]) / dy
        divb = dbx_dx + dby_dy
    else:  # 1D
        divb = (bx[1:] - bx[:-1]) / dx
    
    return divb

def read_athena_binary(filename):
    """
    Read AthenaK binary output file
    
    Returns:
    --------
    data : dict
        Dictionary containing mesh info and variables
    """
    data = {}
    
    try:
        # Check if this is a multi-rank output
        if 'rank_' in filename:
            # Use the rank 0 file to read all ranks
            if 'rank_00000000' in filename:
                filedata = read_all_ranks_binary(filename)
            else:
                # Find rank 0 file
                rank0_file = filename.replace(os.path.basename(filename).split('rank_')[1].split('/')[0], '00000000')
                filedata = read_all_ranks_binary(rank0_file)
        else:
            filedata = read_binary(filename)
        
        # Extract basic info
        data['time'] = filedata['time']
        data['cycle'] = filedata['cycle']
        
        # Extract grid info
        data['x1min'] = filedata['x1min']
        data['x1max'] = filedata['x1max']
        data['x2min'] = filedata['x2min']
        data['x2max'] = filedata['x2max']
        data['x3min'] = filedata['x3min']
        data['x3max'] = filedata['x3max']
        
        # Calculate grid spacing (assuming uniform)
        data['dx'] = (data['x1max'] - data['x1min']) / filedata['Nx1']
        data['dy'] = (data['x2max'] - data['x2min']) / filedata['Nx2'] if filedata['Nx2'] > 1 else 1.0
        data['dz'] = (data['x3max'] - data['x3min']) / filedata['Nx3'] if filedata['Nx3'] > 1 else 1.0
        
        # Create coordinate arrays
        data['x1'] = np.linspace(data['x1min'], data['x1max'], filedata['Nx1'])
        data['x2'] = np.linspace(data['x2min'], data['x2max'], filedata['Nx2'])
        data['x3'] = np.linspace(data['x3min'], data['x3max'], filedata['Nx3'])
        
        # Look for div(B) in the variable list
        var_names = filedata['var_names']
        if 'divb' in var_names or 'mhd_divb' in var_names:
            # Find the index of divb
            divb_idx = None
            for i, name in enumerate(var_names):
                if 'divb' in name.lower():
                    divb_idx = i
                    break
            
            if divb_idx is not None:
                # Extract div(B) data from all mesh blocks and reconstruct full grid
                # For now, just use the first mesh block as a simple test
                if len(filedata['mb_data'][var_names[divb_idx]]) > 0:
                    # Get the first block's data
                    first_block = filedata['mb_data'][var_names[divb_idx]][0]
                    data['divB'] = first_block
                    print(f"  Found div(B) data with shape {first_block.shape}")
        
        # Store raw filedata for additional processing if needed
        data['filedata'] = filedata
        
    except Exception as e:
        print(f"Error reading binary file {filename}: {e}")
        return None
    
    return data

def analyze_divb_statistics(divb, label=""):
    """
    Compute and print div(B) statistics
    """
    divb_abs = np.abs(divb)
    
    print(f"\n{label} div(B) Statistics:")
    print(f"  Mean |div(B)|:     {np.mean(divb_abs):.3e}")
    print(f"  Max |div(B)|:      {np.max(divb_abs):.3e}")
    print(f"  Median |div(B)|:   {np.median(divb_abs):.3e}")
    print(f"  99th percentile:   {np.percentile(divb_abs, 99):.3e}")
    print(f"  99.9th percentile: {np.percentile(divb_abs, 99.9):.3e}")
    
    # Count cells above thresholds
    thresholds = [1e-10, 1e-8, 1e-6, 1e-4]
    total_cells = divb.size
    for thresh in thresholds:
        count = np.sum(divb_abs > thresh)
        percent = 100.0 * count / total_cells
        print(f"  Cells with |div(B)| > {thresh:.0e}: {count} ({percent:.2f}%)")
    
    return divb_abs

def plot_divb_comparison(test_dirs, output_dir):
    """
    Create comparison plots of div(B) for different MPI rank counts
    """
    # Count valid test directories with data
    valid_dirs = {}
    for test_name, test_path in test_dirs.items():
        divb_files = sorted(glob.glob(os.path.join(test_path, "*.divb.*.bin")))
        if not divb_files:
            divb_files = sorted(glob.glob(os.path.join(test_path, "bin", "*.divb.*.bin")))
        if not divb_files:
            divb_files = sorted(glob.glob(os.path.join(test_path, "rank_*/", "*.divb.*.bin")))
        if divb_files:
            valid_dirs[test_name] = test_path
    
    if not valid_dirs:
        print("No valid test directories with div(B) data found for plotting")
        return
    
    n_plots = min(len(valid_dirs), 4)
    fig, axes = plt.subplots(2, n_plots, figsize=(5*n_plots, 10))
    
    if n_plots == 1:
        axes = axes.reshape(2, 1)
    
    max_divb_all = 0
    
    for i, (test_name, test_path) in enumerate(list(valid_dirs.items())[:n_plots]):
        # Find div(B) output files
        divb_files = sorted(glob.glob(os.path.join(test_path, "*.divb.*.bin")))
        if not divb_files:
            divb_files = sorted(glob.glob(os.path.join(test_path, "bin", "*.divb.*.bin")))
        if not divb_files:
            divb_files = sorted(glob.glob(os.path.join(test_path, "rank_*/", "*.divb.*.bin")))
        
        if not divb_files:
            continue
            
        # Use last output file
        data = read_athena_binary(divb_files[-1])
        
        if data is None or 'divB' not in data:
            print(f"div(B) not found in {divb_files[-1]}, skipping...")
            axes[0, i].set_visible(False)
            axes[1, i].set_visible(False)
            continue
        
        divb = data['divB']
        
        # Get dimensions
        if divb.ndim == 3:
            # Take a slice through the middle for 3D
            divb_slice = divb[divb.shape[0]//2, :, :]
            extent = [data['x1'][0], data['x1'][-1], data['x2'][0], data['x2'][-1]]
            xlabel, ylabel = 'x', 'y'
        elif divb.ndim == 2:
            divb_slice = divb
            extent = [data['x1'][0], data['x1'][-1], data['x2'][0], data['x2'][-1]]
            xlabel, ylabel = 'x', 'y'
        else:  # 1D
            axes[0, i].plot(data['x1'], divb)
            axes[0, i].set_xlabel('x')
            axes[0, i].set_ylabel('div(B)')
            axes[0, i].set_title(f'{test_name}\nTime={data["time"]:.2f}')
            
            axes[1, i].semilogy(data['x1'], np.abs(divb))
            axes[1, i].set_xlabel('x')
            axes[1, i].set_ylabel('|div(B)|')
            axes[1, i].set_title('Log scale')
            continue
        
        # Update max for consistent color scale
        max_divb_all = max(max_divb_all, np.max(np.abs(divb_slice)))
        
        # Linear scale plot
        # Use a minimum scale if all values are zero
        vmax = max(max_divb_all, 1e-10)
        im1 = axes[0, i].imshow(divb_slice, extent=extent, origin='lower',
                               cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        axes[0, i].set_xlabel(xlabel)
        axes[0, i].set_ylabel(ylabel)
        axes[0, i].set_title(f'{test_name}\nTime={data["time"]:.2f}')
        plt.colorbar(im1, ax=axes[0, i], label='div(B)')
        
        # Log scale plot
        divb_abs = np.abs(divb_slice)
        divb_abs[divb_abs == 0] = 1e-16  # Avoid log(0)
        
        # Use reasonable bounds for log scale
        vmax_log = max(1e-8, max_divb_all) if max_divb_all > 0 else 1e-8
        im2 = axes[1, i].imshow(divb_abs, extent=extent, origin='lower',
                               cmap='viridis', norm=LogNorm(vmin=1e-12, vmax=vmax_log))
        axes[1, i].set_xlabel(xlabel)
        axes[1, i].set_ylabel(ylabel)
        axes[1, i].set_title('|div(B)| (log scale)')
        plt.colorbar(im2, ax=axes[1, i], label='|div(B)|')
    
    # Hide any unused axes
    for i in range(n_plots, 4):
        if i < axes.shape[1]:
            axes[0, i].set_visible(False)
            axes[1, i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'divb_comparison.png'), dpi=150)
    print(f"\nPlot saved to {os.path.join(output_dir, 'divb_comparison.png')}")

def main(output_dir):
    """
    Main analysis function
    """
    print("=" * 60)
    print("div(B) Analysis for MHD+AMR Tests")
    print("=" * 60)
    
    # Find test directories
    test_dirs = {}
    
    # Check if output_dir itself is a test directory
    if '_np' in os.path.basename(output_dir):
        # Single test directory provided
        test_dirs[os.path.basename(output_dir)] = output_dir
    else:
        # Parent directory provided - look for test subdirectories
        for dirname in os.listdir(output_dir):
            if os.path.isdir(os.path.join(output_dir, dirname)) and '_np' in dirname:
                test_dirs[dirname] = os.path.join(output_dir, dirname)
    
    if not test_dirs:
        print(f"No test directories found in {output_dir}")
        print("Expected directories with format: testname_npN")
        print("\nNote: This script requires binary output files (.bin extension)")
        print("The current test uses binary format (file_type = bin)")
        return
    
    print(f"\nFound {len(test_dirs)} test directories:")
    for name in sorted(test_dirs.keys()):
        print(f"  - {name}")
    
    # Analyze each test
    all_stats = {}
    for test_name in sorted(test_dirs.keys()):
        test_path = test_dirs[test_name]
        print(f"\n{'='*60}")
        print(f"Analyzing: {test_name}")
        print(f"{'='*60}")
        
        # Find output files - look for binary files with .divb. pattern
        divb_files = sorted(glob.glob(os.path.join(test_path, "*.divb.*.bin")))
        if not divb_files:
            # Check in bin subdirectory
            divb_files = sorted(glob.glob(os.path.join(test_path, "bin", "*.divb.*.bin")))
        if not divb_files:
            # Also check in rank subdirectories for multi-rank outputs
            divb_files = sorted(glob.glob(os.path.join(test_path, "rank_*/", "*.divb.*.bin")))
        
        if not divb_files:
            # Check what files exist
            vtk_files = glob.glob(os.path.join(test_path, "vtk/*.vtk"))
            bin_files = glob.glob(os.path.join(test_path, "*.bin"))
            rank_bin_files = glob.glob(os.path.join(test_path, "rank_*/*.bin"))
            
            if vtk_files:
                print(f"  Found VTK files but this script requires binary (.bin) files")
                print(f"  The input file has been updated to use 'file_type = bin'")
                print(f"  Please re-run the test")
            elif bin_files or rank_bin_files:
                all_bins = bin_files + rank_bin_files
                print(f"  Found {len(all_bins)} binary files but none with 'divb' in the name")
                if all_bins:
                    print(f"  Available files: {[os.path.basename(f) for f in all_bins[:3]]}")
                    print(f"  Make sure output3 with variable = mhd_divb is configured") 
            else:
                print(f"  No output files found")
            continue
        
        print(f"  Found {len(divb_files)} output files")
        
        # Analyze first and last outputs
        for idx, fname in enumerate([divb_files[0], divb_files[-1]] if len(divb_files) > 1 else [divb_files[0]]):
            if idx == 0:
                time_label = "Initial"
            else:
                time_label = "Final"
                
            data = read_athena_binary(fname)
            
            if data and 'divB' in data:
                divb = data['divB']
                stats = analyze_divb_statistics(divb, 
                    f"{test_name} - {time_label} (t={data['time']:.2f})")
                all_stats[f"{test_name}_{time_label}"] = stats
            else:
                print(f"  Cannot compute div(B) - data not available in file")
    
    # Compare statistics between different rank counts
    print(f"\n{'='*60}")
    print("Comparison Summary")
    print(f"{'='*60}")
    
    # Group by test type
    test_types = {}
    for key in all_stats.keys():
        base_name = key.split('_np')[0]
        if base_name not in test_types:
            test_types[base_name] = {}
        test_types[base_name][key] = all_stats[key]
    
    for test_type, tests in test_types.items():
        print(f"\n{test_type}:")
        
        # Find single rank baseline
        single_rank = None
        for name in tests.keys():
            if '_np1_' in name:
                single_rank = name
                break
        
        if single_rank:
            baseline_max = np.max(tests[single_rank])
            print(f"  Baseline (1 rank): max|div(B)| = {baseline_max:.3e}")
            
            # Compare others to baseline
            for name in sorted(tests.keys()):
                if name != single_rank:
                    max_divb = np.max(tests[name])
                    ratio = max_divb / baseline_max if baseline_max > 0 else 0
                    nranks = name.split('_np')[1].split('_')[0]
                    print(f"  {nranks} ranks: max|div(B)| = {max_divb:.3e} (ratio = {ratio:.2f})")
    
    # Create comparison plots
    if len(test_dirs) > 0:
        print("\nGenerating comparison plots...")
        plot_divb_comparison(test_dirs, output_dir)
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Check if multiple test directories were provided
        if len(sys.argv) > 2 and all('_np' in arg for arg in sys.argv[1:]):
            # Multiple test directories provided - use parent directory
            output_dir = os.path.dirname(os.path.abspath(sys.argv[1]))
            print(f"Multiple test directories provided, using parent: {output_dir}")
        else:
            output_dir = sys.argv[1]
    else:
        output_dir = "test_divb_output"
    
    if not os.path.exists(output_dir):
        print(f"Error: Directory {output_dir} does not exist")
        print(f"Usage: {sys.argv[0]} [output_directory]")
        print(f"   or: {sys.argv[0]} test_quick_np*")
        sys.exit(1)
    
    main(output_dir)