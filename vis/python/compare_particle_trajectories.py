#!/usr/bin/env python3
"""
Compare particle trajectories before and after restart.

This script reads VTK particle files and compares the positions and properties
of particles from a continuous run vs. a restarted run to verify that restart
preserves particle trajectories correctly.

Usage:
    # Compare all particles in files 2-5
    python compare_particle_trajectories.py --dir-continuous pvtk --dir-restart rst_1/pvtk --fstart 2 --fend 5
    
    # Compare specific particles (tags 100, 200, 300) in files 2-5
    python compare_particle_trajectories.py --dir-continuous pvtk --dir-restart rst_1/pvtk --fstart 2 --fend 5 --ptag 100 200 300
"""

import numpy as np
import struct
import argparse
import sys
from pathlib import Path


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


def compare_particle_states(data1, data2, label1='Run 1', label2='Run 2', 
                           tolerance=1e-6, ptag_filter=None):
    """
    Compare two particle states and report differences.
    
    Args:
        data1: Dictionary from read_vtk_particle_file for first dataset
        data2: Dictionary from read_vtk_particle_file for second dataset
        label1: Label for first dataset
        label2: Label for second dataset
        tolerance: Position difference tolerance
        ptag_filter: List of particle tags to filter (None = compare all particles)
    
    Returns:
        bool: True if states match within tolerance
    """
    print(f"\n{'='*80}")
    print(f"Comparing {label1} vs {label2}")
    if ptag_filter is not None:
        print(f"Filtering for particle tags: {ptag_filter}")
    print(f"{'='*80}")
    
    # Compare metadata
    print(f"\nMetadata comparison:")
    print(f"  {label1}: time={data1['time']:.6e}, cycle={data1['cycle']}, nparticles={data1['nparticles']}")
    print(f"  {label2}: time={data2['time']:.6e}, cycle={data2['cycle']}, nparticles={data2['nparticles']}")
    
    # Only check particle count if not filtering by ptag (will check after filtering)
    if ptag_filter is None and data1['nparticles'] != data2['nparticles']:
        print(f"  ❌ DIFFERENT number of particles!")
        return False
    
    if abs(data1['time'] - data2['time']) > tolerance:
        print(f"  ❌ DIFFERENT times!")
        return False
    
    if ptag_filter is None:
        print(f"  ✓ Metadata matches")
    
    # Sort particles by ptag to ensure same ordering
    if data1['ptag'] is not None and data2['ptag'] is not None:
        idx1 = np.argsort(data1['ptag'])
        idx2 = np.argsort(data2['ptag'])
        
        pos1 = data1['positions'][idx1]
        pos2 = data2['positions'][idx2]
        ptag1 = data1['ptag'][idx1]
        ptag2 = data2['ptag'][idx2]
        gid1 = data1['gid'][idx1] if data1['gid'] is not None else None
        gid2 = data2['gid'][idx2] if data2['gid'] is not None else None
        
        # Filter by ptag if requested
        if ptag_filter is not None:
            mask1 = np.isin(ptag1, ptag_filter)
            mask2 = np.isin(ptag2, ptag_filter)
            
            pos1 = pos1[mask1]
            pos2 = pos2[mask2]
            ptag1 = ptag1[mask1]
            ptag2 = ptag2[mask2]
            gid1 = gid1[mask1] if gid1 is not None else None
            gid2 = gid2[mask2] if gid2 is not None else None
            
            print(f"\n  Filtered to {len(ptag1)} particles in {label1}, {len(ptag2)} particles in {label2}")
            
            if len(ptag1) == 0 or len(ptag2) == 0:
                print(f"  ❌ No particles found with requested tags!")
                return False
        
        # Check if we have the same number of particles after filtering  
        if len(ptag1) != len(ptag2):
            print(f"  ❌ DIFFERENT number of particles!")
            return False
        
        # Check if ptags match
        if not np.array_equal(ptag1, ptag2):
            print(f"  ❌ DIFFERENT particle tags!")
            n_diff = np.sum(ptag1 != ptag2)
            print(f"     {n_diff} particles have different tags")
            return False
        
        print(f"  ✓ Particle tags match ({len(ptag1)} particles compared)")
    else:
        # No tags, compare by position order
        pos1 = data1['positions']
        pos2 = data2['positions']
        gid1 = data1['gid']
        gid2 = data2['gid']
        
        if ptag_filter is not None:
            print(f"  ⚠ Warning: ptag filter requested but no PTAG data in files")
    
    # Compare positions
    pos_diff = np.abs(pos1 - pos2)
    max_diff = np.max(pos_diff)
    mean_diff = np.mean(pos_diff)
    rms_diff = np.sqrt(np.mean(pos_diff**2))
    
    print(f"\nPosition differences:")
    print(f"  Max difference:  {max_diff:.6e}")
    print(f"  Mean difference: {mean_diff:.6e}")
    print(f"  RMS difference:  {rms_diff:.6e}")
    
    if max_diff > tolerance:
        print(f"  ❌ Position differences exceed tolerance {tolerance:.6e}")
        # Find worst offenders
        worst_idx = np.argmax(np.max(pos_diff, axis=1))
        print(f"\n  Worst particle (index {worst_idx}):")
        if data1['ptag'] is not None and data2['ptag'] is not None:
            print(f"    PTAG: {ptag1[worst_idx]}")
        print(f"    {label1}: {pos1[worst_idx]}")
        print(f"    {label2}: {pos2[worst_idx]}")
        print(f"    Difference: {pos_diff[worst_idx]}")
        return False
    else:
        print(f"  ✓ Positions match within tolerance")
    
    # Compare GIDs if available
    if gid1 is not None and gid2 is not None:
        gid_match = np.array_equal(gid1, gid2)
        if gid_match:
            print(f"  ✓ Global IDs match")
        else:
            n_diff = np.sum(gid1 != gid2)
            print(f"  ⚠ {n_diff} particles have different GIDs")
            print(f"    (This is expected if particles moved between MeshBlocks)")
    
    print(f"\n✅ States match within tolerance!")
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Compare particle trajectories before and after restart'
    )
    parser.add_argument('--dir-continuous', type=str, required=True,
                        help='Directory containing VTK particle files from continuous run')
    parser.add_argument('--dir-restart', type=str, required=True,
                        help='Directory containing VTK particle files from restarted run')
    parser.add_argument('--fstart', type=int, required=True,
                        help='First file number to compare (e.g., 1 for KH.prtcl_all.00001.part.vtk)')
    parser.add_argument('--fend', type=int, required=True,
                        help='Last file number to compare (e.g., 5 for KH.prtcl_all.00005.part.vtk)')
    parser.add_argument('--basename', type=str, default='KH.prtcl_all',
                        help='Base name of particle files (default: KH.prtcl_all)')
    parser.add_argument('--tolerance', type=float, default=1e-6,
                        help='Position difference tolerance (default: 1e-6)')
    parser.add_argument('--ptag', type=int, nargs='+',
                        help='Specific particle tag(s) to compare (optional, compare all if not specified)')
    
    args = parser.parse_args()
    
    dir_continuous = Path(args.dir_continuous)
    dir_restart = Path(args.dir_restart)
    
    if not dir_continuous.exists():
        print(f"Error: Directory {dir_continuous} does not exist")
        sys.exit(1)
    
    if not dir_restart.exists():
        print(f"Error: Directory {dir_restart} does not exist")
        sys.exit(1)
    
    fstart = args.fstart
    fend = args.fend
    
    if fstart > fend:
        print(f"Error: fstart ({fstart}) must be <= fend ({fend})")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print(f"Particle Trajectory Comparison")
    print(f"{'='*80}")
    print(f"Continuous run directory: {dir_continuous}")
    print(f"Restarted run directory:  {dir_restart}")
    print(f"Comparing files {fstart} to {fend}")
    print(f"Position tolerance: {args.tolerance:.6e}")
    if args.ptag is not None:
        print(f"Filtering for particle tags: {args.ptag}")
    
    all_match = True
    
    # Compare timesteps in the specified range
    for step_num in range(fstart, fend + 1):
        
        # Construct filenames
        file_continuous = dir_continuous / f"{args.basename}.{step_num:05d}.part.vtk"
        file_restart = dir_restart / f"{args.basename}.{step_num:05d}.part.vtk"
        
        if not file_continuous.exists():
            print(f"\n⚠ Warning: File {file_continuous} does not exist, stopping comparison")
            break
        
        if not file_restart.exists():
            print(f"\n⚠ Warning: File {file_restart} does not exist, stopping comparison")
            break
        
        print(f"\n{'─'*80}")
        print(f"Step {step_num}")
        print(f"{'─'*80}")
        
        try:
            data_continuous = read_vtk_particle_file(file_continuous)
            data_restart = read_vtk_particle_file(file_restart)
            
            matches = compare_particle_states(
                data_continuous, 
                data_restart,
                label1='Continuous',
                label2='Restarted',
                tolerance=args.tolerance,
                ptag_filter=args.ptag
            )
            
            if not matches:
                all_match = False
        
        except Exception as e:
            print(f"❌ Error comparing files: {e}")
            import traceback
            traceback.print_exc()
            all_match = False
    
    print(f"\n{'='*80}")
    if all_match:
        print("✅ SUCCESS: All compared timesteps match!")
        print("   Restart preserves particle trajectories correctly.")
    else:
        print("❌ FAILURE: Some timesteps do not match!")
        print("   Restart may not be preserving particle trajectories correctly.")
    print(f"{'='*80}\n")
    
    return 0 if all_match else 1


if __name__ == '__main__':
    sys.exit(main())
