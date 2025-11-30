#!/usr/bin/env python3
"""
Read AthenaK HDF5 particle output files (.h5)

This module provides functions to read particle data stored in HDF5 format.
The HDF5 files are typically created by converting binary particle files using
prtcl_bin_to_hdf5.py

Usage:
    import read_prtcl_hdf5 as rph
    data = rph.read_particle_hdf5('pbin/TRML.00001.h5')
    
    # Access data
    print(f"Time: {data['time']}, Cycle: {data['ncycle']}")
    print(f"Number of particles: {data['nparticles']}")
    print(f"Particle positions: x={data['x']}, y={data['y']}, z={data['z']}")
"""

import numpy as np
import h5py
import sys


def read_particle_hdf5(filename):
    """
    Read AthenaK HDF5 particle output file
    
    Parameters:
    -----------
    filename : str
        Path to the HDF5 particle file (.h5)
        
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
    
    with h5py.File(filename, 'r') as f:
        # Read metadata
        metadata = f['Metadata']
        data['time'] = metadata.attrs['time']
        data['dt'] = metadata.attrs['dt']
        data['ncycle'] = metadata.attrs['ncycle']
        data['nparticles'] = metadata.attrs['nparticles']
        data['nrdata'] = metadata.attrs['nrdata']
        data['nidata'] = metadata.attrs['nidata']
        data['ngriddata'] = metadata.attrs['ngriddata']
        
        # Read variable names
        var_names_bytes = metadata['var_names'][:]
        data['var_names'] = [vn.decode('utf-8') for vn in var_names_bytes]
        
        # Read particle real data
        if 'ParticleRealData' in f:
            rdata_grp = f['ParticleRealData']
            for key in rdata_grp.keys():
                data[key] = rdata_grp[key][:]
        
        # Read particle integer data
        if 'ParticleIntData' in f:
            idata_grp = f['ParticleIntData']
            for key in idata_grp.keys():
                data[key] = idata_grp[key][:]
        
        # Read grid data
        if 'GridData' in f:
            griddata_grp = f['GridData']
            for key in griddata_grp.keys():
                data[key] = griddata_grp[key][:]
    
    return data


def list_hdf5_contents(filename):
    """
    Print the structure and contents of an HDF5 particle file
    
    Parameters:
    -----------
    filename : str
        Path to the HDF5 particle file (.h5)
    """
    
    print(f"HDF5 file: {filename}")
    print("="*70)
    
    with h5py.File(filename, 'r') as f:
        # Print metadata
        if 'Metadata' in f:
            print("\nMetadata:")
            metadata = f['Metadata']
            for attr_name, attr_value in metadata.attrs.items():
                print(f"  {attr_name}: {attr_value}")
            
            if 'var_names' in metadata:
                var_names = [vn.decode('utf-8') for vn in metadata['var_names'][:]]
                print(f"  var_names: {var_names}")
        
        # Print particle real data
        if 'ParticleRealData' in f:
            print("\nParticle Real Data:")
            rdata_grp = f['ParticleRealData']
            for key in rdata_grp.keys():
                dset = rdata_grp[key]
                print(f"  {key}: shape={dset.shape}, dtype={dset.dtype}")
        
        # Print particle integer data
        if 'ParticleIntData' in f:
            print("\nParticle Integer Data:")
            idata_grp = f['ParticleIntData']
            for key in idata_grp.keys():
                dset = idata_grp[key]
                print(f"  {key}: shape={dset.shape}, dtype={dset.dtype}")
        
        # Print grid data
        if 'GridData' in f:
            print("\nGrid Data at Particle Locations:")
            griddata_grp = f['GridData']
            for key in griddata_grp.keys():
                dset = griddata_grp[key]
                print(f"  {key}: shape={dset.shape}, dtype={dset.dtype}")
    
    print("="*70)


def print_particle_info(data):
    """
    Print summary information about the particle data
    
    Parameters:
    -----------
    data : dict
        Dictionary returned by read_particle_hdf5()
    """
    print("="*70)
    print("AthenaK HDF5 Particle Data Summary")
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
        if varname in data:
            print(f"  {varname}: [{data[varname].min():.6e}, {data[varname].max():.6e}]")
    
    print("\nParticle metadata:")
    print(f"  GID range: [{data['gid'].min()}, {data['gid'].max()}]")
    print(f"  Tag range: [{data['tag'].min()}, {data['tag'].max()}]")
    print("="*70)


def read_particle_subset(filename, indices=None, variables=None):
    """
    Read a subset of particles or variables from HDF5 file
    
    This is useful for large files where you only need specific particles
    or specific variables.
    
    Parameters:
    -----------
    filename : str
        Path to the HDF5 particle file (.h5)
    indices : array-like, optional
        Indices of particles to read. If None, read all particles.
    variables : list of str, optional
        List of variable names to read. If None, read all variables.
        
    Returns:
    --------
    data : dict
        Dictionary containing requested data
    """
    
    data = {}
    
    with h5py.File(filename, 'r') as f:
        # Always read metadata
        metadata = f['Metadata']
        data['time'] = metadata.attrs['time']
        data['dt'] = metadata.attrs['dt']
        data['ncycle'] = metadata.attrs['ncycle']
        data['nparticles'] = metadata.attrs['nparticles']
        data['nrdata'] = metadata.attrs['nrdata']
        data['nidata'] = metadata.attrs['nidata']
        data['ngriddata'] = metadata.attrs['ngriddata']
        
        var_names_bytes = metadata['var_names'][:]
        data['var_names'] = [vn.decode('utf-8') for vn in var_names_bytes]
        
        # Determine which variables to read
        if variables is None:
            read_all_vars = True
        else:
            read_all_vars = False
            variables = set(variables)
        
        # Read particle real data
        if 'ParticleRealData' in f:
            rdata_grp = f['ParticleRealData']
            for key in rdata_grp.keys():
                if read_all_vars or key in variables:
                    if indices is None:
                        data[key] = rdata_grp[key][:]
                    else:
                        data[key] = rdata_grp[key][indices]
        
        # Read particle integer data
        if 'ParticleIntData' in f:
            idata_grp = f['ParticleIntData']
            for key in idata_grp.keys():
                if read_all_vars or key in variables:
                    if indices is None:
                        data[key] = idata_grp[key][:]
                    else:
                        data[key] = idata_grp[key][indices]
        
        # Read grid data
        if 'GridData' in f:
            griddata_grp = f['GridData']
            for key in griddata_grp.keys():
                if read_all_vars or key in variables:
                    if indices is None:
                        data[key] = griddata_grp[key][:]
                    else:
                        data[key] = griddata_grp[key][indices]
    
    # Update nparticles if reading subset
    if indices is not None:
        data['nparticles'] = len(indices)
    
    return data


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python read_prtcl_hdf5.py <filename.h5> [--list]")
        print("\nExamples:")
        print("  python read_prtcl_hdf5.py pbin/TRML.00001.h5")
        print("  python read_prtcl_hdf5.py pbin/TRML.00001.h5 --list")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    if '--list' in sys.argv:
        list_hdf5_contents(filename)
    else:
        try:
            data = read_particle_hdf5(filename)
            print_particle_info(data)
        except Exception as e:
            print(f"Error reading file: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
