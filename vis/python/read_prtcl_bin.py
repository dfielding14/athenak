#!/usr/bin/env python3
"""
Read AthenaK binary particle output files (.prtclbin)

Binary file format:
  1. Header (variable size):
     - int64 magic_number (= 43)
     - int64 total_particles
     - int64 nrdata (number of real data fields per particle)
     - int64 nidata (number of integer data fields per particle)
     - int64 ngriddata (number of grid data fields per particle)
     - double time
     - double dt
     - double ncycle
     - char[16] var_names[ngriddata] (variable names, 16 chars each, null-padded)
  2. Particle real data (nrdata arrays of nprtcl doubles each)
  3. Particle integer data (nidata arrays of nprtcl int32 each)
  4. Grid data at particle locations (ngriddata arrays of nprtcl doubles each)

Usage:
    import read_prtcl_bin as rpb
    data = rpb.read_particle_binary('pbin/TRML.00001.prtclbin')
    
    # Access data
    print(f"Time: {data['time']}, Cycle: {data['ncycle']}")
    print(f"Number of particles: {data['nparticles']}")
    print(f"Particle positions: x={data['x']}, y={data['y']}, z={data['z']}")
    print(f"Particle velocities: vx={data['vx']}, vy={data['vy']}, vz={data['vz']}")
    print(f"Particle tags: {data['tag']}, gids: {data['gid']}")
    print(f"Grid data variable names: {data['var_names']}")
    # Grid data keys depend on var_names, e.g.:
    print(f"Density at particles: {data['dens']}")
    print(f"Temperature at particles: {data['temp']}")
"""

import numpy as np
import struct
import sys

def read_particle_binary(filename):
    """
    Read AthenaK binary particle output file
    
    Parameters:
    -----------
    filename : str
        Path to the binary particle file (.prtclbin)
        
    Returns:
    --------
    data : dict
        Dictionary containing:
        - 'time': simulation time
        - 'dt': timestep
        - 'ncycle': cycle number
        - 'nparticles': total number of particles
        - 'nrdata': number of real data fields
        - 'nidata': number of integer data fields
        - 'ngriddata': number of grid data fields
        - 'var_names': list of grid variable names
        - 'x', 'y', 'z': particle positions
        - 'vx', 'vy', 'vz': particle velocities (for Lagrangian tracers)
        - 'gid': MeshBlock global IDs
        - 'tag': particle tags
        - <var_name>: grid data at particle locations (names from var_names list)
    """
    
    data = {}
    
    with open(filename, 'rb') as f:
        # Read header (variable size)
        # 5 int64_t values (40 bytes)
        header_int = np.fromfile(f, dtype=np.int64, count=5)
        magic_number = header_int[0]
        nparticles = header_int[1]
        nrdata = header_int[2]
        nidata = header_int[3]
        ngriddata = header_int[4]
        
        # 3 double values (24 bytes)
        header_real = np.fromfile(f, dtype=np.float64, count=3)
        time = header_real[0]
        dt = header_real[1]
        ncycle = int(header_real[2])
        
        # Verify magic number
        if magic_number != 43:
            raise ValueError(f"Invalid magic number {magic_number}, expected 43")
        
        # Read variable names (16 bytes each)
        var_names = []
        for i in range(ngriddata):
            varname_bytes = f.read(16)
            # Decode and strip null bytes
            varname = varname_bytes.decode('utf-8').rstrip('\x00')
            var_names.append(varname)
        
        # Store header info
        data['time'] = time
        data['dt'] = dt
        data['ncycle'] = ncycle
        data['nparticles'] = nparticles
        data['nrdata'] = nrdata
        data['nidata'] = nidata
        data['ngriddata'] = ngriddata
        data['var_names'] = var_names
        
        # Read particle real data
        # For Lagrangian tracers: x, y, z, vx, vy, vz
        rdata_names = ['x', 'y', 'z', 'vx', 'vy', 'vz']
        for i in range(nrdata):
            rdata = np.fromfile(f, dtype=np.float64, count=nparticles)
            if i < len(rdata_names):
                data[rdata_names[i]] = rdata
            else:
                data[f'rdata_{i}'] = rdata
        
        # Read particle integer data
        # Typically: gid, tag, lastmove, lastlevel
        idata_names = ['gid', 'tag', 'lastmove', 'lastlevel']
        for i in range(nidata):
            idata = np.fromfile(f, dtype=np.int32, count=nparticles)
            if i < len(idata_names):
                data[idata_names[i]] = idata
            else:
                data[f'idata_{i}'] = idata
        
        # Read grid data at particle locations
        # Use variable names from file header
        for i in range(ngriddata):
            griddata = np.fromfile(f, dtype=np.float64, count=nparticles)
            data[var_names[i]] = griddata
    
    return data


def print_particle_info(data):
    """
    Print summary information about the particle data
    
    Parameters:
    -----------
    data : dict
        Dictionary returned by read_particle_binary()
    """
    print("="*70)
    print("AthenaK Binary Particle Data Summary")
    print("="*70)
    print(f"Time:            {data['time']:.6e}")
    print(f"Timestep (dt):   {data['dt']:.6e}")
    print(f"Cycle number:    {data['ncycle']}")
    print(f"Total particles: {data['nparticles']}")
    print(f"Real data fields: {data['nrdata']}")
    print(f"Int data fields:  {data['nidata']}")
    print(f"Grid data fields: {data['ngriddata']}")
    print(f"Variable names:  {data['var_names']}")
    print("="*70)
    
    print("\nParticle position ranges:")
    print(f"  x: [{data['x'].min():.6e}, {data['x'].max():.6e}]")
    print(f"  y: [{data['y'].min():.6e}, {data['y'].max():.6e}]")
    print(f"  z: [{data['z'].min():.6e}, {data['z'].max():.6e}]")
    
    if 'vx' in data:
        print("\nParticle velocity ranges:")
        print(f"  vx: [{data['vx'].min():.6e}, {data['vx'].max():.6e}]")
        print(f"  vy: [{data['vy'].min():.6e}, {data['vy'].max():.6e}]")
        print(f"  vz: [{data['vz'].min():.6e}, {data['vz'].max():.6e}]")
    
    print("\nGrid data at particle locations:")
    for varname in data['var_names']:
        print(f"  {varname}: [{data[varname].min():.6e}, {data[varname].max():.6e}]")
    
    print("\nParticle metadata:")
    print(f"  GID range: [{data['gid'].min()}, {data['gid'].max()}]")
    print(f"  Tag range: [{data['tag'].min()}, {data['tag'].max()}]")
    print("="*70)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python read_prtcl_bin.py <filename.prtclbin>")
        print("\nExample:")
        print("  python read_prtcl_bin.py pbin/TRML.00001.prtclbin")
        sys.exit(1)
    
    filename = sys.argv[1]
    print(f"Reading particle data from: {filename}\n")
    
    try:
        data = read_particle_binary(filename)
        print_particle_info(data)
        
        # Optional: save to npz for further analysis
        output_npz = filename.replace('.prtclbin', '.npz')
        np.savez(output_npz, **data)
        print(f"\nData also saved to: {output_npz}")
        print("Load with: data = np.load('{output_npz}')")
        
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
