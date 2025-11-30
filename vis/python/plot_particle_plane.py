#!/usr/bin/env python3
"""
Plot particle data projected onto a 2D plane.

This script reads VTK particle files and creates 2D projections onto xy, xz, or yz planes.
Particles are colored by time using cmasher colormaps, allowing visualization of 
particle evolution and comparison between different simulation runs.

Usage:
    # Plot xy plane for single simulation
    python plot_particle_plane.py --dir pvtk --plane xy --fstart 1 --fend 10 --output particles_xy.png
    
    # Compare multiple simulations with different colormaps
    python plot_particle_plane.py --dir pvtk rst_1/pvtk rst_2/pvtk --plane xz --fstart 1 --fend 10 \
           --labels "No Restart" "Restart 1" "Restart 2" --cmaps cmr.ember cmr.ocean cmr.forest \
           --output comparison.png
    
    # Plot specific particle tags
    python plot_particle_plane.py --dir pvtk --plane yz --fstart 1 --fend 10 --ptag 100 200 300 \
           --output tagged_particles.png
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import argparse
import sys
from pathlib import Path

# Try to import cmasher, fall back to matplotlib colormaps if not available
try:
    import cmasher as cmr
    HAS_CMASHER = True
except ImportError:
    HAS_CMASHER = False
    print("Warning: cmasher not found. Using matplotlib colormaps instead.")
    print("Install with: pip install cmasher")


def read_vtk_particle_file(filename):
    """
    Read a legacy VTK particle file and extract particle data.
    
    Returns:
        dict: Dictionary containing:
            - 'time': simulation time
            - 'cycle': cycle number
            - 'nparticles': number of particles
            - 'positions': array of (x,y,z) positions, shape (N,3)
            - 'gid': global ID array, shape (N,)
            - 'ptag': particle tag array, shape (N,)
    """
    with open(filename, 'rb') as f:
        # Read ASCII header
        header_lines = []
        while True:
            line = f.readline().decode('ascii')
            header_lines.append(line.strip())
            if 'POINTS' in line:
                # Extract number of particles
                parts = line.split()
                nparticles = int(parts[1])
                break
            if len(header_lines) > 100:  # Safety check
                raise ValueError(f"Could not find POINTS line in {filename}")
        
        # Extract time and cycle from header comment
        time = None
        cycle = None
        for line in header_lines:
            if 'time=' in line:
                # Parse: "# AthenaK particle data at time= 0.086349  nranks= 4  cycle=281  variables=..."
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == 'time=' and i+1 < len(parts):
                        time = float(parts[i+1])
                    if part.startswith('cycle='):
                        cycle = int(part.split('=')[1])
        
        # Read binary particle positions (3*nparticles floats in big-endian)
        position_data = f.read(4 * 3 * nparticles)
        positions = np.frombuffer(position_data, dtype='>f4').reshape((nparticles, 3))
        
        # Read POINT_DATA header
        while True:
            line = f.readline().decode('ascii')
            if 'POINT_DATA' in line:
                break
            if not line:
                raise ValueError(f"Could not find POINT_DATA in {filename}")
        
        # Read scalar data (gid, ptag, etc.)
        gid = None
        ptag = None
        
        while True:
            line = f.readline()
            if not line:
                break
            line = line.decode('ascii').strip()
            
            if 'SCALARS gid' in line:
                # Skip LOOKUP_TABLE line
                f.readline()
                # Read gid data
                gid_data = f.read(4 * nparticles)
                gid = np.frombuffer(gid_data, dtype='>f4').astype(np.int32)
            
            elif 'SCALARS ptag' in line:
                # Skip LOOKUP_TABLE line
                f.readline()
                # Read ptag data
                ptag_data = f.read(4 * nparticles)
                ptag = np.frombuffer(ptag_data, dtype='>f4').astype(np.int32)
    
    return {
        'time': time,
        'cycle': cycle,
        'nparticles': nparticles,
        'positions': positions,
        'gid': gid,
        'ptag': ptag
    }


def get_colormap(cmap_name):
    """
    Get a colormap, trying cmasher first, then falling back to matplotlib.
    
    Args:
        cmap_name: Name of the colormap (e.g., 'cmr.ember', 'viridis')
    
    Returns:
        matplotlib colormap object
    """
    if HAS_CMASHER and cmap_name.startswith('cmr.'):
        # Use cmasher colormap
        cmap_short_name = cmap_name[4:]  # Remove 'cmr.' prefix
        try:
            return cmr.__dict__[cmap_short_name]
        except KeyError:
            print(f"Warning: cmasher colormap '{cmap_short_name}' not found, using 'viridis'")
            return plt.cm.viridis
    else:
        # Use matplotlib colormap
        try:
            return plt.cm.get_cmap(cmap_name)
        except:
            print(f"Warning: colormap '{cmap_name}' not found, using 'viridis'")
            return plt.cm.viridis


def plot_particle_plane(directories, plane='xy', fstart=1, fend=10, basename='KH.prtcl_all',
                        labels=None, cmaps=None, ptag_filter=None, output='particles.png',
                        figsize=(12, 10), alpha=0.5, marker_size=1, xlim=None, ylim=None,
                        verbose=False):
    """
    Plot particle positions projected onto a 2D plane.
    
    Args:
        directories: List of directory paths containing VTK particle files
        plane: Projection plane ('xy', 'xz', or 'yz')
        fstart: First file number to read
        fend: Last file number to read
        basename: Base name of particle files
        labels: List of labels for each directory (for legend)
        cmaps: List of colormap names for each directory
        ptag_filter: List of particle tags to plot (None = plot all)
        output: Output filename for the plot
        figsize: Figure size (width, height) in inches
        alpha: Transparency of particle markers
        marker_size: Size of particle markers
        xlim: Tuple of (xmin, xmax) for x-axis limits (None = auto)
        ylim: Tuple of (ymin, ymax) for y-axis limits (None = auto)
    verbose: If True, print detailed particle position information for tagged particles
    """
    # Convert single directory to list
    if isinstance(directories, (str, Path)):
        directories = [directories]
    
    # Set default labels if not provided
    if labels is None:
        labels = [f"Run {i+1}" for i in range(len(directories))]
    
    # Set default colormaps if not provided
    if cmaps is None:
        if HAS_CMASHER:
            default_cmaps = ['cmr.ember', 'cmr.ocean', 'cmr.forest', 'cmr.lavender', 
                           'cmr.cosmic', 'cmr.dusk', 'cmr.arctic', 'cmr.pepper']
        else:
            default_cmaps = ['viridis', 'plasma', 'inferno', 'magma', 
                           'cividis', 'twilight', 'turbo', 'jet']
        cmaps = [default_cmaps[i % len(default_cmaps)] for i in range(len(directories))]
    
    # Validate plane selection
    plane = plane.lower()
    if plane not in ['xy', 'xz', 'yz']:
        raise ValueError(f"Invalid plane '{plane}'. Must be 'xy', 'xz', or 'yz'")
    
    # Determine axis indices
    axis_map = {
        'xy': (0, 1, 'X', 'Y'),
        'xz': (0, 2, 'X', 'Z'),
        'yz': (1, 2, 'Y', 'Z')
    }
    idx1, idx2, label1, label2 = axis_map[plane]
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Track global min/max for auto-scaling
    global_x_min, global_x_max = None, None
    global_y_min, global_y_max = None, None
    
    # Process each directory
    for dir_idx, (directory, label, cmap_name) in enumerate(zip(directories, labels, cmaps)):
        directory = Path(directory)
        
        if not directory.exists():
            print(f"Warning: Directory {directory} does not exist, skipping")
            continue
        
        print(f"\nProcessing {label} ({directory})...")
        
        # Get colormap
        cmap = get_colormap(cmap_name)
        
        # Collect all particle positions and times
        all_positions = []
        all_times = []
        
        for step_num in range(fstart, fend + 1):
            filename = directory / f"{basename}.{step_num:05d}.part.vtk"
            
            if not filename.exists():
                print(f"  Warning: {filename} not found, stopping at step {step_num-1}")
                break
            
            try:
                data = read_vtk_particle_file(filename)
                
                # Filter by ptag if requested
                if ptag_filter is not None and data['ptag'] is not None:
                    mask = np.isin(data['ptag'], ptag_filter)
                    positions = data['positions'][mask]
                    filtered_ptags = data['ptag'][mask]
                    n_filtered = np.sum(mask)
                    print(f"  Step {step_num}: time={data['time']:.6e}, "
                          f"{n_filtered}/{data['nparticles']} particles (filtered)")
                    
                    # Verbose output: print individual particle positions
                    if verbose and n_filtered > 0:
                        print(f"    Tagged particle positions in {label}:")
                        # Sort by ptag for consistent ordering
                        sort_idx = np.argsort(filtered_ptags)
                        for idx in sort_idx:
                            ptag = filtered_ptags[idx]
                            pos = positions[idx]
                            print(f"      PTAG {ptag:6d}: x={pos[0]:18.15f}, y={pos[1]:18.15f}, z={pos[2]:18.15f}")
                    
                    # (analysis removed) track per-step positions if requested via separate tool
                else:
                    positions = data['positions']
                    print(f"  Step {step_num}: time={data['time']:.6e}, "
                          f"{data['nparticles']} particles")
                
                # Store positions and times
                all_positions.append(positions)
                all_times.append(np.full(len(positions), data['time']))
                
            except Exception as e:
                print(f"  Error reading {filename}: {e}")
                continue
        
        if len(all_positions) == 0:
            print(f"  Warning: No valid data found for {label}")
            continue
        
        # analysis-summary removed from plotting script (use analyze_particle_movement.py)
        
        # Concatenate all positions and times
        all_positions = np.vstack(all_positions)
        all_times = np.concatenate(all_times)
        
        # Extract coordinates for the selected plane
        x = all_positions[:, idx1]
        y = all_positions[:, idx2]
        
        # Update global min/max for auto-scaling
        if global_x_min is None:
            global_x_min = np.min(x)
            global_x_max = np.max(x)
            global_y_min = np.min(y)
            global_y_max = np.max(y)
        else:
            global_x_min = min(global_x_min, np.min(x))
            global_x_max = max(global_x_max, np.max(x))
            global_y_min = min(global_y_min, np.min(y))
            global_y_max = max(global_y_max, np.max(y))
        
        # Normalize times for colormap
        time_min = np.min(all_times)
        time_max = np.max(all_times)
        if time_max > time_min:
            norm = Normalize(vmin=time_min, vmax=time_max)
            colors = cmap(norm(all_times))
        else:
            colors = cmap(0.5)  # Single color if all times are the same
        
        # Plot particles
        scatter = ax.scatter(x, y, c=all_times, cmap=cmap, s=marker_size, 
                           alpha=alpha, label=label, edgecolors='none')
        
        # Add colorbar for the first dataset (or create separate colorbars)
        if dir_idx == 0 or len(directories) == 1:
            cbar = plt.colorbar(scatter, ax=ax, label='Time')
            cbar.ax.tick_params(labelsize=10)
    
    # Set labels and title
    ax.set_xlabel(label1, fontsize=14)
    ax.set_ylabel(label2, fontsize=14)
    ax.set_title(f'Particle Distribution - {plane.upper()} Plane', fontsize=16, fontweight='bold')
    
    # Set axis limits
    # Use provided limits if given, otherwise use data min/max
    if xlim is not None:
        ax.set_xlim(xlim)
        print(f"\nUsing user-specified x-axis limits: [{xlim[0]:.4f}, {xlim[1]:.4f}]")
    elif global_x_min is not None:
        # Add small padding (2% of range) for better visualization
        x_range = global_x_max - global_x_min
        padding = 0.02 * x_range if x_range > 0 else 0.1
        ax.set_xlim(global_x_min - padding, global_x_max + padding)
        print(f"\nAuto-detected {label1}-axis limits: [{global_x_min:.4f}, {global_x_max:.4f}]")
        print(f"  (with 2% padding: [{global_x_min - padding:.4f}, {global_x_max + padding:.4f}])")
    
    if ylim is not None:
        ax.set_ylim(ylim)
        print(f"Using user-specified y-axis limits: [{ylim[0]:.4f}, {ylim[1]:.4f}]")
    elif global_y_min is not None:
        # Add small padding (2% of range) for better visualization
        y_range = global_y_max - global_y_min
        padding = 0.02 * y_range if y_range > 0 else 0.1
        ax.set_ylim(global_y_min - padding, global_y_max + padding)
        print(f"Auto-detected {label2}-axis limits: [{global_y_min:.4f}, {global_y_max:.4f}]")
        print(f"  (with 2% padding: [{global_y_min - padding:.4f}, {global_y_max + padding:.4f}])")
    
    # Add legend if multiple datasets
    if len(directories) > 1:
        ax.legend(loc='upper right', fontsize=10, framealpha=0.8)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Equal aspect ratio for physical coordinates
    ax.set_aspect('equal')
    
    # Tight layout
    plt.tight_layout()
    
    # Save figure
    print(f"\nSaving plot to {output}...")
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Plot saved successfully!")
    
    # Show plot
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Plot particle positions projected onto a 2D plane',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single simulation, xy plane
  python plot_particle_plane.py --dir pvtk --plane xy --fstart 1 --fend 10
  
  # Compare two simulations with different colormaps
  python plot_particle_plane.py --dir pvtk rst_1/pvtk --plane xz --fstart 1 --fend 10 \\
         --labels "Continuous" "Restarted" --cmaps cmr.ember cmr.ocean
  
  # Plot specific particles
  python plot_particle_plane.py --dir pvtk --plane yz --fstart 1 --fend 10 \\
         --ptag 100 200 300 --output tagged_particles.png
  
Available cmasher colormaps (if installed):
  cmr.ember, cmr.ocean, cmr.forest, cmr.lavender, cmr.cosmic, cmr.dusk,
  cmr.arctic, cmr.pepper, cmr.flamingo, cmr.voltage, cmr.freeze, cmr.heat
  
Fallback matplotlib colormaps:
  viridis, plasma, inferno, magma, cividis, twilight, turbo, jet
        """
    )
    
    parser.add_argument('--dir', type=str, nargs='+', required=True,
                        help='Directory(ies) containing VTK particle files')
    parser.add_argument('--plane', type=str, choices=['xy', 'xz', 'yz'], default='xy',
                        help='Projection plane (default: xy)')
    parser.add_argument('--fstart', type=int, required=True,
                        help='First file number to read')
    parser.add_argument('--fend', type=int, required=True,
                        help='Last file number to read')
    parser.add_argument('--basename', type=str, default='KH.prtcl_all',
                        help='Base name of particle files (default: KH.prtcl_all)')
    parser.add_argument('--labels', type=str, nargs='+',
                        help='Labels for each directory (for legend)')
    parser.add_argument('--cmaps', type=str, nargs='+',
                        help='Colormap names for each directory (e.g., cmr.ember viridis)')
    parser.add_argument('--ptag', type=int, nargs='+',
                        help='Specific particle tag(s) to plot (optional, plot all if not specified)')
    parser.add_argument('--output', type=str, default='particles.png',
                        help='Output filename (default: particles.png)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[12, 10],
                        help='Figure size in inches: width height (default: 12 10)')
    parser.add_argument('--alpha', type=float, default=0.5,
                        help='Transparency of particle markers (default: 0.5)')
    parser.add_argument('--marker-size', type=float, default=1.0,
                        help='Size of particle markers (default: 1.0)')
    parser.add_argument('--xlim', type=float, nargs=2,
                        help='X-axis limits: xmin xmax')
    parser.add_argument('--ylim', type=float, nargs=2,
                        help='Y-axis limits: ymin ymax')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print detailed particle position information for tagged particles')
    # analysis flags moved to scripts/analyze_particle_movement.py
    args = parser.parse_args()
    
    # Validate inputs
    if args.fstart > args.fend:
        print(f"Error: fstart ({args.fstart}) must be <= fend ({args.fend})")
        sys.exit(1)
    
    if args.labels is not None and len(args.labels) != len(args.dir):
        print(f"Error: Number of labels ({len(args.labels)}) must match number of directories ({len(args.dir)})")
        sys.exit(1)
    
    if args.cmaps is not None and len(args.cmaps) != len(args.dir):
        print(f"Error: Number of colormaps ({len(args.cmaps)}) must match number of directories ({len(args.dir)})")
        sys.exit(1)
    
    # Create plot
    plot_particle_plane(
        directories=args.dir,
        plane=args.plane,
        fstart=args.fstart,
        fend=args.fend,
        basename=args.basename,
        labels=args.labels,
        cmaps=args.cmaps,
        ptag_filter=args.ptag,
        output=args.output,
        figsize=tuple(args.figsize),
        alpha=args.alpha,
        marker_size=args.marker_size,
        xlim=tuple(args.xlim) if args.xlim else None,
        ylim=tuple(args.ylim) if args.ylim else None,
        verbose=args.verbose
    )
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
