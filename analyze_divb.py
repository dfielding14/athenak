#!/usr/bin/env python3
"""
analyze_divb.py
Analyze div(B) from AthenaK MHD+AMR simulation outputs
Compares div(B) errors between single-rank and multi-rank runs
"""

import sys
import os
import numpy as np
import h5py
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm

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

def read_athena_hdf5(filename):
    """
    Read AthenaK HDF5 output file
    
    Returns:
    --------
    data : dict
        Dictionary containing mesh info and variables
    """
    data = {}
    
    with h5py.File(filename, 'r') as f:
        # Read mesh information
        data['time'] = f.attrs['Time']
        data['cycle'] = f.attrs['Cycle']
        
        # Read coordinate information
        data['x1'] = f['x1v'][:]
        data['x2'] = f['x2v'][:]
        data['x3'] = f['x3v'][:]
        
        data['x1f'] = f['x1f'][:]
        data['x2f'] = f['x2f'][:]
        data['x3f'] = f['x3f'][:]
        
        # Calculate grid spacing
        data['dx'] = data['x1'][1] - data['x1'][0] if len(data['x1']) > 1 else 1.0
        data['dy'] = data['x2'][1] - data['x2'][0] if len(data['x2']) > 1 else 1.0
        data['dz'] = data['x3'][1] - data['x3'][0] if len(data['x3']) > 1 else 1.0
        
        # Read magnetic field if present
        if 'B' in f:
            # Cell-centered B (for some outputs)
            data['Bcc'] = f['B'][:]
        
        # Read div(B) if directly output
        if 'divB' in f:
            data['divB'] = f['divB'][:]
        elif 'mhd_divb' in f:
            data['divB'] = f['mhd_divb'][:]
        
        # Read face-centered fields if available
        if 'Bx1f' in f:
            data['Bx1f'] = f['Bx1f'][:]
        if 'Bx2f' in f:
            data['Bx2f'] = f['Bx2f'][:]
        if 'Bx3f' in f:
            data['Bx3f'] = f['Bx3f'][:]
            
        # Read MeshBlock information if available
        if 'MeshBlockIndex' in f:
            data['mb_index'] = f['MeshBlockIndex'][:]
            
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
    fig, axes = plt.subplots(2, len(test_dirs), figsize=(5*len(test_dirs), 10))
    
    if len(test_dirs) == 1:
        axes = axes.reshape(2, 1)
    
    max_divb_all = 0
    
    for i, (test_name, test_path) in enumerate(test_dirs.items()):
        # Find div(B) output files
        divb_files = sorted(glob.glob(os.path.join(test_path, "*divb*.athdf")))
        
        if not divb_files:
            print(f"No div(B) files found in {test_path}")
            continue
            
        # Use last output file
        data = read_athena_hdf5(divb_files[-1])
        
        if 'divB' in data:
            divb = data['divB']
        else:
            print(f"div(B) not found in {divb_files[-1]}, skipping...")
            continue
        
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
        im1 = axes[0, i].imshow(divb_slice, extent=extent, origin='lower',
                               cmap='RdBu_r', vmin=-max_divb_all, vmax=max_divb_all)
        axes[0, i].set_xlabel(xlabel)
        axes[0, i].set_ylabel(ylabel)
        axes[0, i].set_title(f'{test_name}\nTime={data["time"]:.2f}')
        plt.colorbar(im1, ax=axes[0, i], label='div(B)')
        
        # Log scale plot
        divb_abs = np.abs(divb_slice)
        divb_abs[divb_abs == 0] = 1e-16  # Avoid log(0)
        
        im2 = axes[1, i].imshow(divb_abs, extent=extent, origin='lower',
                               cmap='viridis', norm=LogNorm(vmin=1e-12, vmax=max(1e-8, max_divb_all)))
        axes[1, i].set_xlabel(xlabel)
        axes[1, i].set_ylabel(ylabel)
        axes[1, i].set_title('|div(B)| (log scale)')
        plt.colorbar(im2, ax=axes[1, i], label='|div(B)|')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'divb_comparison.png'), dpi=150)
    print(f"\nPlot saved to {os.path.join(output_dir, 'divb_comparison.png')}")
    plt.show()

def main(output_dir):
    """
    Main analysis function
    """
    print("=" * 60)
    print("div(B) Analysis for MHD+AMR Tests")
    print("=" * 60)
    
    # Find test directories
    test_dirs = {}
    for dirname in os.listdir(output_dir):
        if os.path.isdir(os.path.join(output_dir, dirname)) and '_np' in dirname:
            test_dirs[dirname] = os.path.join(output_dir, dirname)
    
    if not test_dirs:
        print(f"No test directories found in {output_dir}")
        print("Expected directories with format: testname_npN")
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
        
        # Find output files
        divb_files = sorted(glob.glob(os.path.join(test_path, "*divb*.athdf")))
        
        if not divb_files:
            print(f"  No div(B) output files found")
            continue
        
        print(f"  Found {len(divb_files)} output files")
        
        # Analyze first and last outputs
        for idx, fname in enumerate([divb_files[0], divb_files[-1]]):
            if idx == 0:
                time_label = "Initial"
            else:
                time_label = "Final"
                
            data = read_athena_hdf5(fname)
            
            if 'divB' in data:
                divb = data['divB']
                stats = analyze_divb_statistics(divb, 
                    f"{test_name} - {time_label} (t={data['time']:.2f})")
                all_stats[f"{test_name}_{time_label}"] = stats
            else:
                # Try to calculate from face-centered fields
                if 'Bx1f' in data and 'Bx2f' in data:
                    print(f"  Calculating div(B) from face-centered fields...")
                    divb = calculate_divb_from_fields(
                        data['Bx1f'], data['Bx2f'], 
                        data.get('Bx3f', np.zeros_like(data['Bx1f'])),
                        data['dx'], data['dy'], data['dz']
                    )
                    stats = analyze_divb_statistics(divb,
                        f"{test_name} - {time_label} (t={data['time']:.2f})")
                    all_stats[f"{test_name}_{time_label}"] = stats
                else:
                    print(f"  Cannot compute div(B) - fields not available")
    
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
        output_dir = sys.argv[1]
    else:
        output_dir = "test_divb_output"
    
    if not os.path.exists(output_dir):
        print(f"Error: Directory {output_dir} does not exist")
        print(f"Usage: {sys.argv[0]} [output_directory]")
        sys.exit(1)
    
    main(output_dir)