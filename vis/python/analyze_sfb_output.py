#!/usr/bin/env python3
"""
Analyze output from SFB turbulence simulations
This script reads AthenaK output files and computes various diagnostics
for spherical Fourier-Bessel turbulence simulations.
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
import glob
import os
from scipy import stats
from matplotlib.animation import FuncAnimation
from matplotlib.colors import SymLogNorm

def read_athena_hdf5(filename):
    """
    Read AthenaK HDF5 output file
    """
    data = {}
    with h5py.File(filename, 'r') as f:
        # Read metadata
        data['time'] = f.attrs['Time']
        data['cycle'] = f.attrs['Cycle']
        
        # Read mesh information
        data['x1v'] = f['x1v'][:]
        data['x2v'] = f['x2v'][:]
        data['x3v'] = f['x3v'][:]
        
        # Read primitive variables
        # Assumes hydro output with density, velocity, pressure
        data['rho'] = f['prim'][:,:,:,0]
        data['vx'] = f['prim'][:,:,:,1]
        data['vy'] = f['prim'][:,:,:,2]
        data['vz'] = f['prim'][:,:,:,3]
        data['press'] = f['prim'][:,:,:,4]
        
    return data

def read_history_file(filename):
    """
    Read AthenaK history file
    """
    # Read the history file
    data = np.loadtxt(filename, skiprows=2)
    
    # Read header to get column names
    with open(filename, 'r') as f:
        lines = f.readlines()
        header = lines[1].strip().split()[1:]  # Skip the # character
    
    # Create dictionary
    history = {}
    for i, col in enumerate(header):
        history[col] = data[:, i]
    
    return history

def compute_spherical_averages(data, r0_turb):
    """
    Compute spherically averaged quantities
    """
    # Create coordinate grids
    x1v, x2v, x3v = np.meshgrid(data['x1v'], data['x2v'], data['x3v'], indexing='ij')
    
    # Compute radius
    r = np.sqrt(x1v**2 + x2v**2 + x3v**2)
    
    # Create radial bins
    r_bins = np.linspace(0, r0_turb*1.2, 50)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    
    # Initialize arrays for averages
    rho_avg = np.zeros_like(r_centers)
    v_rms = np.zeros_like(r_centers)
    v_rad = np.zeros_like(r_centers)
    mach_avg = np.zeros_like(r_centers)
    
    # Compute velocity magnitude
    v_mag = np.sqrt(data['vx']**2 + data['vy']**2 + data['vz']**2)
    
    # Compute radial velocity
    v_r = (data['vx']*x1v + data['vy']*x2v + data['vz']*x3v) / (r + 1e-10)
    
    # Compute sound speed (assuming gamma = 5/3)
    gamma = 5.0/3.0
    cs = np.sqrt(gamma * data['press'] / data['rho'])
    mach = v_mag / cs
    
    # Bin by radius
    for i in range(len(r_centers)):
        mask = (r >= r_bins[i]) & (r < r_bins[i+1])
        if np.any(mask):
            rho_avg[i] = np.mean(data['rho'][mask])
            v_rms[i] = np.sqrt(np.mean(v_mag[mask]**2))
            v_rad[i] = np.mean(v_r[mask])
            mach_avg[i] = np.mean(mach[mask])
    
    return {
        'r': r_centers,
        'rho': rho_avg,
        'v_rms': v_rms,
        'v_rad': v_rad,
        'mach': mach_avg,
        'r_3d': r,
        'v_mag_3d': v_mag,
        'mach_3d': mach
    }

def compute_power_spectrum(data):
    """
    Compute power spectrum of velocity field
    """
    # Get grid spacing
    dx = data['x1v'][1] - data['x1v'][0]
    
    # Compute FFT of each velocity component
    vx_fft = np.fft.fftn(data['vx'])
    vy_fft = np.fft.fftn(data['vy'])
    vz_fft = np.fft.fftn(data['vz'])
    
    # Compute power
    power = np.abs(vx_fft)**2 + np.abs(vy_fft)**2 + np.abs(vz_fft)**2
    
    # Get wavenumbers
    nx = len(data['x1v'])
    kx = np.fft.fftfreq(nx, d=dx) * 2 * np.pi
    
    # Create k grid
    kx_grid, ky_grid, kz_grid = np.meshgrid(kx, kx, kx, indexing='ij')
    k_mag = np.sqrt(kx_grid**2 + ky_grid**2 + kz_grid**2)
    
    # Bin by k magnitude
    k_bins = np.logspace(0, 2.5, 40)
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
    power_spectrum = np.zeros(len(k_centers))
    
    for i in range(len(k_centers)):
        mask = (k_mag >= k_bins[i]) & (k_mag < k_bins[i+1])
        if np.any(mask):
            power_spectrum[i] = np.mean(power[mask])
    
    return k_centers, power_spectrum

def check_divergence(data):
    """
    Check divergence of velocity field
    """
    dx = data['x1v'][1] - data['x1v'][0]
    dy = data['x2v'][1] - data['x2v'][0]
    dz = data['x3v'][1] - data['x3v'][0]
    
    # Compute divergence using central differences
    dvx_dx = np.gradient(data['vx'], dx, axis=0)
    dvy_dy = np.gradient(data['vy'], dy, axis=1)
    dvz_dz = np.gradient(data['vz'], dz, axis=2)
    
    div = dvx_dx + dvy_dy + dvz_dz
    
    return div

def plot_slice(data, r0_turb, output_file=None):
    """
    Plot slices through the simulation
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Get middle slice
    nz = len(data['x3v'])
    iz = nz // 2
    
    # Create coordinate grids for plotting
    x1v_2d, x2v_2d = np.meshgrid(data['x1v'], data['x2v'], indexing='ij')
    
    # Density
    im1 = axes[0,0].pcolormesh(x1v_2d, x2v_2d, data['rho'][:,:,iz].T, 
                               cmap='viridis', shading='auto')
    axes[0,0].set_title('Density')
    axes[0,0].set_xlabel('x')
    axes[0,0].set_ylabel('y')
    plt.colorbar(im1, ax=axes[0,0])
    
    # Velocity magnitude
    v_mag = np.sqrt(data['vx']**2 + data['vy']**2 + data['vz']**2)
    im2 = axes[0,1].pcolormesh(x1v_2d, x2v_2d, v_mag[:,:,iz].T, 
                               cmap='plasma', shading='auto')
    axes[0,1].set_title('Velocity Magnitude')
    axes[0,1].set_xlabel('x')
    axes[0,1].set_ylabel('y')
    plt.colorbar(im2, ax=axes[0,1])
    
    # Velocity vectors
    skip = max(len(data['x1v']) // 20, 1)
    axes[0,2].quiver(x1v_2d[::skip,::skip], x2v_2d[::skip,::skip],
                     data['vx'][::skip,::skip,iz].T, data['vy'][::skip,::skip,iz].T,
                     alpha=0.7)
    axes[0,2].set_title('Velocity Vectors')
    axes[0,2].set_xlabel('x')
    axes[0,2].set_ylabel('y')
    axes[0,2].set_aspect('equal')
    
    # Vorticity (z-component)
    dx = data['x1v'][1] - data['x1v'][0]
    dy = data['x2v'][1] - data['x2v'][0]
    dvy_dx = np.gradient(data['vy'][:,:,iz], dx, axis=0)
    dvx_dy = np.gradient(data['vx'][:,:,iz], dy, axis=1)
    vorticity = dvy_dx - dvx_dy
    
    im4 = axes[1,0].pcolormesh(x1v_2d, x2v_2d, vorticity.T, 
                               cmap='RdBu_r', shading='auto',
                               norm=SymLogNorm(linthresh=0.1))
    axes[1,0].set_title('Vorticity (z-component)')
    axes[1,0].set_xlabel('x')
    axes[1,0].set_ylabel('y')
    plt.colorbar(im4, ax=axes[1,0])
    
    # Divergence
    div = check_divergence(data)
    im5 = axes[1,1].pcolormesh(x1v_2d, x2v_2d, div[:,:,iz].T, 
                               cmap='RdBu_r', shading='auto',
                               norm=SymLogNorm(linthresh=1e-10))
    axes[1,1].set_title('Divergence')
    axes[1,1].set_xlabel('x')
    axes[1,1].set_ylabel('y')
    plt.colorbar(im5, ax=axes[1,1])
    
    # Mach number
    gamma = 5.0/3.0
    cs = np.sqrt(gamma * data['press'] / data['rho'])
    mach = v_mag / cs
    im6 = axes[1,2].pcolormesh(x1v_2d, x2v_2d, mach[:,:,iz].T, 
                               cmap='hot', shading='auto')
    axes[1,2].set_title('Mach Number')
    axes[1,2].set_xlabel('x')
    axes[1,2].set_ylabel('y')
    plt.colorbar(im6, ax=axes[1,2])
    
    # Add circles showing r0_turb
    for ax in axes.flat:
        circle = plt.Circle((0, 0), r0_turb, fill=False, color='white', 
                          linewidth=2, linestyle='--', alpha=0.7)
        ax.add_patch(circle)
        ax.set_aspect('equal')
    
    plt.suptitle(f'SFB Turbulence at t = {data["time"]:.3f}, cycle = {data["cycle"]}')
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    else:
        plt.show()

def plot_radial_profiles(spherical_avg, r0_turb, output_file=None):
    """
    Plot spherically averaged radial profiles
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    r = spherical_avg['r']
    
    # Density profile
    axes[0,0].plot(r, spherical_avg['rho'], 'b-', linewidth=2)
    axes[0,0].axvline(r0_turb, color='r', linestyle='--', alpha=0.7, label='r0_turb')
    axes[0,0].set_xlabel('r')
    axes[0,0].set_ylabel('⟨ρ⟩')
    axes[0,0].set_title('Density Profile')
    axes[0,0].grid(True, alpha=0.3)
    axes[0,0].legend()
    
    # RMS velocity profile
    axes[0,1].plot(r, spherical_avg['v_rms'], 'g-', linewidth=2)
    axes[0,1].axvline(r0_turb, color='r', linestyle='--', alpha=0.7, label='r0_turb')
    axes[0,1].set_xlabel('r')
    axes[0,1].set_ylabel('v_rms')
    axes[0,1].set_title('RMS Velocity Profile')
    axes[0,1].grid(True, alpha=0.3)
    axes[0,1].legend()
    
    # Radial velocity profile
    axes[1,0].plot(r, spherical_avg['v_rad'], 'r-', linewidth=2)
    axes[1,0].axvline(r0_turb, color='r', linestyle='--', alpha=0.7, label='r0_turb')
    axes[1,0].axhline(0, color='k', linestyle='-', alpha=0.3)
    axes[1,0].set_xlabel('r')
    axes[1,0].set_ylabel('⟨v_r⟩')
    axes[1,0].set_title('Mean Radial Velocity')
    axes[1,0].grid(True, alpha=0.3)
    axes[1,0].legend()
    
    # Mach number profile
    axes[1,1].plot(r, spherical_avg['mach'], 'm-', linewidth=2)
    axes[1,1].axvline(r0_turb, color='r', linestyle='--', alpha=0.7, label='r0_turb')
    axes[1,1].axhline(1, color='k', linestyle=':', alpha=0.5, label='M=1')
    axes[1,1].set_xlabel('r')
    axes[1,1].set_ylabel('⟨M⟩')
    axes[1,1].set_title('Mach Number Profile')
    axes[1,1].grid(True, alpha=0.3)
    axes[1,1].legend()
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    else:
        plt.show()

def plot_power_spectrum(k, power, expected_slope=-5/3, output_file=None):
    """
    Plot power spectrum with expected slope
    """
    plt.figure(figsize=(8, 6))
    
    # Plot power spectrum
    plt.loglog(k, power, 'b-', linewidth=2, label='Measured')
    
    # Plot expected slope
    k_ref = k[len(k)//3]
    power_ref = power[len(power)//3]
    expected = power_ref * (k/k_ref)**expected_slope
    plt.loglog(k, expected, 'r--', linewidth=2, 
               label=f'k^{{{expected_slope:.2f}}}')
    
    plt.xlabel('k')
    plt.ylabel('E(k)')
    plt.title('Velocity Power Spectrum')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    else:
        plt.show()

def analyze_time_series(history, output_dir=None):
    """
    Analyze time series from history file
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    time = history['time']
    
    # Kinetic energy
    if 'Ek/V' in history:
        axes[0,0].plot(time, history['Ek/V'], 'b-', linewidth=2)
        axes[0,0].set_xlabel('Time')
        axes[0,0].set_ylabel('Ek/V')
        axes[0,0].set_title('Volume-averaged Kinetic Energy')
        axes[0,0].grid(True, alpha=0.3)
    
    # Mach number
    if 'Mach' in history:
        axes[0,1].plot(time, history['Mach'], 'g-', linewidth=2)
        axes[0,1].set_xlabel('Time')
        axes[0,1].set_ylabel('Mach')
        axes[0,1].set_title('Mass-weighted Mach Number')
        axes[0,1].grid(True, alpha=0.3)
    
    # Density fluctuations
    if 'drho_rms/rho' in history:
        axes[1,0].plot(time, history['drho_rms/rho'], 'r-', linewidth=2)
        axes[1,0].set_xlabel('Time')
        axes[1,0].set_ylabel('δρ/ρ')
        axes[1,0].set_title('RMS Density Fluctuations')
        axes[1,0].grid(True, alpha=0.3)
    
    # Time step
    if 'dt' in history:
        axes[1,1].plot(time, history['dt'], 'm-', linewidth=2)
        axes[1,1].set_xlabel('Time')
        axes[1,1].set_ylabel('dt')
        axes[1,1].set_title('Time Step')
        axes[1,1].grid(True, alpha=0.3)
        axes[1,1].set_yscale('log')
    
    plt.tight_layout()
    
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'time_series.png'), 
                   dpi=150, bbox_inches='tight')
    else:
        plt.show()

def print_statistics(data, spherical_avg, div, r0_turb):
    """
    Print key statistics
    """
    print("\n" + "="*60)
    print(f"SFB Turbulence Statistics at t = {data['time']:.3f}")
    print("="*60)
    
    # Get mask for turbulent region
    r = spherical_avg['r_3d']
    mask = r < r0_turb
    
    # Velocity statistics
    v_mag = spherical_avg['v_mag_3d']
    print("\nVelocity Field (inside r0_turb):")
    print(f"  Mean |v|: {np.mean(v_mag[mask]):.4f}")
    print(f"  RMS velocity: {np.sqrt(np.mean(v_mag[mask]**2)):.4f}")
    print(f"  Max |v|: {np.max(v_mag[mask]):.4f}")
    
    # Momentum conservation
    print("\nMomentum Conservation:")
    print(f"  <vx>: {np.mean(data['vx'][mask]):.6e}")
    print(f"  <vy>: {np.mean(data['vy'][mask]):.6e}")
    print(f"  <vz>: {np.mean(data['vz'][mask]):.6e}")
    print(f"  |<v>|: {np.sqrt(np.mean(data['vx'][mask])**2 + np.mean(data['vy'][mask])**2 + np.mean(data['vz'][mask])**2):.6e}")
    
    # Divergence statistics
    print("\nDivergence Statistics (inside r0_turb):")
    print(f"  Max |div(v)|: {np.max(np.abs(div[mask])):.6e}")
    print(f"  Mean |div(v)|: {np.mean(np.abs(div[mask])):.6e}")
    print(f"  RMS div(v): {np.sqrt(np.mean(div[mask]**2)):.6e}")
    
    # Mach number
    mach = spherical_avg['mach_3d']
    print("\nMach Number (inside r0_turb):")
    print(f"  Mean Mach: {np.mean(mach[mask]):.4f}")
    print(f"  Max Mach: {np.max(mach[mask]):.4f}")
    
    # Density statistics
    rho = data['rho']
    rho_mean = np.mean(rho[mask])
    print("\nDensity (inside r0_turb):")
    print(f"  Mean ρ: {rho_mean:.4f}")
    print(f"  δρ/ρ (RMS): {np.std(rho[mask])/rho_mean:.4f}")
    
    print("="*60 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Analyze SFB turbulence output')
    parser.add_argument('input', help='Input file or pattern (e.g., "output/*.athdf")')
    parser.add_argument('--r0_turb', type=float, default=0.4,
                        help='Turbulence region radius (default: 0.4)')
    parser.add_argument('--history', help='History file (e.g., "sfb_turb.hst")')
    parser.add_argument('--output-dir', help='Output directory for plots')
    parser.add_argument('--all-frames', action='store_true',
                        help='Analyze all frames matching pattern')
    parser.add_argument('--power-spectrum', action='store_true',
                        help='Compute and plot power spectrum')
    parser.add_argument('--expo', type=float, default=5/3,
                        help='Expected power spectrum slope (default: 5/3)')
    
    args = parser.parse_args()
    
    # Create output directory if specified
    if args.output_dir and not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Handle file patterns
    if '*' in args.input:
        files = sorted(glob.glob(args.input))
        if not files:
            print(f"No files matching pattern: {args.input}")
            return
    else:
        files = [args.input]
    
    # Analyze history file if provided
    if args.history and os.path.exists(args.history):
        print(f"Analyzing history file: {args.history}")
        history = read_history_file(args.history)
        analyze_time_series(history, args.output_dir)
    
    # Analyze data files
    for i, filename in enumerate(files):
        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            continue
        
        print(f"\nAnalyzing: {filename}")
        
        # Read data
        data = read_athena_hdf5(filename)
        
        # Compute spherical averages
        spherical_avg = compute_spherical_averages(data, args.r0_turb)
        
        # Check divergence
        div = check_divergence(data)
        
        # Print statistics
        print_statistics(data, spherical_avg, div, args.r0_turb)
        
        # Generate output filename base
        if args.output_dir:
            base_name = os.path.splitext(os.path.basename(filename))[0]
            slice_file = os.path.join(args.output_dir, f'{base_name}_slice.png')
            radial_file = os.path.join(args.output_dir, f'{base_name}_radial.png')
        else:
            slice_file = None
            radial_file = None
        
        # Plot slices
        if not args.all_frames or i == 0:  # Plot first frame or all
            plot_slice(data, args.r0_turb, slice_file)
            plot_radial_profiles(spherical_avg, args.r0_turb, radial_file)
        
        # Compute power spectrum for last frame
        if args.power_spectrum and i == len(files) - 1:
            print("Computing power spectrum...")
            k, power = compute_power_spectrum(data)
            
            if args.output_dir:
                spectrum_file = os.path.join(args.output_dir, 'power_spectrum.png')
            else:
                spectrum_file = None
            
            plot_power_spectrum(k, power, -args.expo, spectrum_file)

if __name__ == '__main__':
    main()