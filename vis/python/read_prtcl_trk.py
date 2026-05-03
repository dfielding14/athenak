#!/usr/bin/env python3
"""
Read AthenaK tracked particle files (.trk)

The .trk format appends data at each output time:
  - Text header line with metadata (time, nranks, cycle, ntracked_prtcls)
  - Binary data: nvalues floats per tracked particle. New files include:
    x, y, z, vx, vy, vz, rho, press, temp, eint, scalar0, gid, level, active
    ordered by particle tag (0 to ntracked_prtcls-1)

This format allows tracking specific particles over time by their tag.

Usage:
    import read_prtcl_trk as rpt
    data = rpt.read_tracked_particles('trk/TRML.trk')
    
    # Access data - list of snapshots, each a dict
    for i, snapshot in enumerate(data):
        print(f"Snapshot {i}: t={snapshot['time']}, cycle={snapshot['cycle']}")
        print(f"  Particle 0 position: {snapshot['positions'][0]}")
        print(f"  Particle 0 velocity: {snapshot['velocities'][0]}")
"""

import numpy as np
import sys
import re


def read_tracked_particles(filename):
    """
    Read AthenaK tracked particle file (.trk)
    
    Parameters:
    -----------
    filename : str
        Path to the .trk file
        
    Returns:
    --------
    snapshots : list of dicts
        List of snapshots, each containing:
        - 'time': simulation time
        - 'nranks': number of MPI ranks
        - 'cycle': cycle number
        - 'ntracked': number of tracked particles
        - 'positions': ndarray of shape (ntracked, 3) with (x,y,z)
        - 'velocities': ndarray of shape (ntracked, 3) with (vx,vy,vz)
        - 'fields': dict of all variables present in the file
    """
    
    snapshots = []
    
    with open(filename, 'rb') as f:
        while True:
            # Try to read a header line
            pos = f.tell()
            line = f.readline()
            
            if not line:
                break  # End of file
                
            try:
                line_str = line.decode('ascii').strip()
            except UnicodeDecodeError:
                # Hit binary data, rewind and skip
                f.seek(pos)
                # Skip some bytes and try again
                f.read(1024)
                continue
            
            # Check if this is a header line
            if '# AthenaK tracked particle data' in line_str:
                # Parse header
                # Format: # AthenaK tracked particle data at time= XXX  nranks= YYY  cycle=ZZZ  ntracked_prtcls=NNN
                match = re.search(r'time=\s*([\d.eE+-]+)', line_str)
                time = float(match.group(1)) if match else None
                
                match = re.search(r'nranks=\s*(\d+)', line_str)
                nranks = int(match.group(1)) if match else None
                
                match = re.search(r'cycle=\s*(\d+)', line_str)
                cycle = int(match.group(1)) if match else None
                
                match = re.search(r'ntracked_prtcls=\s*(\d+)', line_str)
                ntracked = int(match.group(1)) if match else None

                match = re.search(r'nvalues=\s*(\d+)', line_str)
                nvalues = int(match.group(1)) if match else 6

                match = re.search(r'variables=\s*(.*)$', line_str)
                if match:
                    variables = match.group(1).split()
                else:
                    variables = ['x', 'y', 'z', 'vx', 'vy', 'vz']

                if len(variables) != nvalues:
                    variables = [f'value{i}' for i in range(nvalues)]
                
                if ntracked is None:
                    print(f"Warning: Could not parse ntracked_prtcls from header: {line_str}")
                    continue
                
                data = np.fromfile(f, dtype=np.float32, count=nvalues*ntracked)
                
                if len(data) < nvalues*ntracked:
                    print(f"Warning: Expected {nvalues*ntracked} values but got {len(data)}")
                    break
                
                data = data.reshape((ntracked, nvalues))
                fields = {name: data[:, i] for i, name in enumerate(variables)}
                
                snapshot = {
                    'time': time,
                    'nranks': nranks,
                    'cycle': cycle,
                    'ntracked': ntracked,
                    'nvalues': nvalues,
                    'variables': variables,
                    'fields': fields,
                    'positions': data[:, 0:3],
                    'velocities': data[:, 3:6]
                }
                
                snapshots.append(snapshot)
    
    return snapshots


def print_tracked_summary(snapshots):
    """
    Print summary of tracked particle data
    
    Parameters:
    -----------
    snapshots : list
        List of snapshot dicts from read_tracked_particles()
    """
    print("="*70)
    print("AthenaK Tracked Particle Data Summary")
    print("="*70)
    print(f"Total snapshots: {len(snapshots)}")
    
    if len(snapshots) == 0:
        print("No data found in file")
        return
    
    first = snapshots[0]
    last = snapshots[-1]
    
    print(f"\nFirst snapshot:")
    print(f"  Time:   {first['time']:.6e}")
    print(f"  Cycle:  {first['cycle']}")
    print(f"  Tracked particles: {first['ntracked']}")
    
    print(f"\nLast snapshot:")
    print(f"  Time:   {last['time']:.6e}")
    print(f"  Cycle:  {last['cycle']}")
    print(f"  Tracked particles: {last['ntracked']}")
    
    print(f"\nTime coverage: {first['time']:.6e} to {last['time']:.6e}")
    print(f"Variables: {' '.join(first.get('variables', []))}")
    
    # Show example particle trajectory
    if first['ntracked'] > 0:
        print(f"\nParticle 0 trajectory:")
        print(f"  Initial position: {first['positions'][0]}")
        print(f"  Final position:   {last['positions'][0]}")
        displacement = last['positions'][0] - first['positions'][0]
        print(f"  Net displacement: {displacement}")
        print(f"  Distance traveled: {np.linalg.norm(displacement):.6e}")
    
    print("="*70)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python read_prtcl_trk.py <filename.trk>")
        print("\nExample:")
        print("  python read_prtcl_trk.py trk/TRML.trk")
        sys.exit(1)
    
    filename = sys.argv[1]
    print(f"Reading tracked particle data from: {filename}\n")
    
    try:
        snapshots = read_tracked_particles(filename)
        print_tracked_summary(snapshots)
        
        # Optional: save to npz
        if len(sys.argv) > 2 and sys.argv[2] == '--save-npz':
            output_npz = filename.replace('.trk', '_tracked.npz')
            # Combine all snapshots into arrays
            times = np.array([s['time'] for s in snapshots])
            cycles = np.array([s['cycle'] for s in snapshots])
            
            # Stack positions and velocities: shape (nsnapshots, ntracked, 3)
            ntracked = snapshots[0]['ntracked']
            positions = np.array([s['positions'] for s in snapshots])
            velocities = np.array([s['velocities'] for s in snapshots])
            field_arrays = {
                name: np.array([s['fields'][name] for s in snapshots])
                for name in snapshots[0].get('variables', [])
            }
            
            np.savez(output_npz,
                     times=times,
                     cycles=cycles,
                     ntracked=ntracked,
                     positions=positions,
                     velocities=velocities,
                     **field_arrays)
            print(f"\nData saved to: {output_npz}")
        
    except Exception as e:
        print(f"Error reading file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
