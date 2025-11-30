#!/usr/bin/env python3
"""
Read binary particle restart files (.prtclrst) produced by AthenaK

Binary file format:
  Header: 2 x int64_t
    - Magic number (42) for format identification
    - Total number of particles
  Data: 6 x float64 arrays, each of length nparticles
    - gid[nparticles]      : MeshBlock global ID
    - tag[nparticles]      : Particle tag
    - plastmove[nparticles]: Last move status (-1 = frozen/deleted, >=0 = active)
    - x[nparticles]        : X position
    - y[nparticles]        : Y position  
    - z[nparticles]        : Z position
"""

import numpy as np
import sys
import os


def read_particle_restart(filename):
    """
    Read a binary particle restart file.
    
    Parameters:
    -----------
    filename : str
        Path to the .prtclrst file
        
    Returns:
    --------
    dict : Dictionary containing:
        'magic' : int - Magic number (should be 42)
        'nparticles' : int - Total number of particles
        'gid' : ndarray - MeshBlock global IDs
        'tag' : ndarray - Particle tags
        'plastmove' : ndarray - Last move status (-1=frozen/deleted, >=0=active)
        'x' : ndarray - X positions
        'y' : ndarray - Y positions
        'z' : ndarray - Z positions
    """
    with open(filename, 'rb') as f:
        # Read header
        header = np.fromfile(f, dtype=np.int64, count=2)
        magic, nparticles = header
        
        if magic != 42:
            print(f"Warning: Magic number is {magic}, expected 42")
        
        # Read particle data
        gid = np.fromfile(f, dtype=np.float64, count=nparticles)
        tag = np.fromfile(f, dtype=np.float64, count=nparticles)
        plastmove = np.fromfile(f, dtype=np.float64, count=nparticles)
        x = np.fromfile(f, dtype=np.float64, count=nparticles)
        y = np.fromfile(f, dtype=np.float64, count=nparticles)
        z = np.fromfile(f, dtype=np.float64, count=nparticles)
        
        # Verify we read all the data
        if len(gid) != nparticles or len(tag) != nparticles or \
           len(plastmove) != nparticles or \
           len(x) != nparticles or len(y) != nparticles or len(z) != nparticles:
            raise ValueError(f"Data mismatch: expected {nparticles} particles")
    
    return {
        'magic': magic,
        'nparticles': nparticles,
        'gid': gid,
        'tag': tag,
        'plastmove': plastmove,
        'x': x,
        'y': y,
        'z': z
    }


def print_particle_summary(data):
    """Print summary information about the particle data."""
    print(f"\n=== Particle Restart File Summary ===")
    print(f"Magic number:      {data['magic']}")
    print(f"Total particles:   {data['nparticles']}")
    print(f"\nGID range:         [{data['gid'].min():.0f}, {data['gid'].max():.0f}]")
    print(f"TAG range:         [{data['tag'].min():.0f}, {data['tag'].max():.0f}]")
    print(f"PLASTMOVE range:   [{data['plastmove'].min():.0f}, {data['plastmove'].max():.0f}]")
    n_frozen = np.sum(data['plastmove'] == -1)
    n_deleted = np.sum(data['plastmove'] == -2)
    n_active = np.sum(data['plastmove'] >= 0)
    print(f"  Active particles:  {n_active}")
    print(f"  Frozen particles:  {n_frozen}")
    print(f"  Deleted particles: {n_deleted}")
    print(f"X range:           [{data['x'].min():.6e}, {data['x'].max():.6e}]")
    print(f"Y range:           [{data['y'].min():.6e}, {data['y'].max():.6e}]")
    print(f"Z range:           [{data['z'].min():.6e}, {data['z'].max():.6e}]")
    
    print(f"\nFirst 5 particles:")
    print(f"  GID:       {data['gid'][:5]}")
    print(f"  TAG:       {data['tag'][:5]}")
    print(f"  PLASTMOVE: {data['plastmove'][:5]}")
    print(f"  X:         {data['x'][:5]}")
    print(f"  Y:         {data['y'][:5]}")
    print(f"  Z:         {data['z'][:5]}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python read_particle_restart.py <filename.prtclrst>")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found")
        sys.exit(1)
    
    # Read the particle data
    data = read_particle_restart(filename)
    
    # Print summary
    print_particle_summary(data)
    
    # Optionally save to a different format
    if len(sys.argv) > 2 and sys.argv[2] == '--save-npz':
        outfile = filename.replace('.prtclrst', '.npz')
        np.savez(outfile, 
                 gid=data['gid'], 
                 tag=data['tag'],
                 plastmove=data['plastmove'],
                 x=data['x'], 
                 y=data['y'], 
                 z=data['z'])
        print(f"\nSaved data to {outfile}")
