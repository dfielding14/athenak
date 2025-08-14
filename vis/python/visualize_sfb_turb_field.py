#!/usr/bin/env python3
"""
Visualize SFB (Spherical Fourier-Bessel) turbulence driving field with AMR
Shows radial structure and angular patterns characteristic of SFB basis
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from bin_convert_new import read_binary_as_athdf

def visualize_sfb_turb_field(filename, output_name):
    """Create visualization of SFB turbulence driving field"""
    
    # Read binary data
    print(f"Reading {filename}...")
    data = read_binary_as_athdf(filename)
    
    # Extract force components
    force_x = data['force1']
    force_y = data['force2']
    force_z = data['force3']
    
    # Get dimensions
    nz, ny, nx = force_x.shape
    print(f"Data dimensions: {nx} x {ny} x {nz}")
    
    # Calculate magnitude of force
    force_mag = np.sqrt(force_x**2 + force_y**2 + force_z**2)
    
    # Create figure with multiple panels showing SFB characteristics
    fig = plt.figure(figsize=(18, 14))
    
    # Get spatial coordinates
    x = np.linspace(-1.0, 1.0, nx)
    y = np.linspace(-1.0, 1.0, ny)
    z = np.linspace(-1.0, 1.0, nz)
    
    # Panel 1: Force magnitude at z=0 (xy plane)
    ax1 = plt.subplot(3, 3, 1)
    k_slice = nz // 2
    fmag_xy = force_mag[k_slice, :, :]
    im1 = ax1.imshow(fmag_xy, origin='lower', cmap='viridis', 
                     extent=[-1, 1, -1, 1], interpolation='bilinear')
    ax1.set_title('Force Magnitude (z=0)', fontsize=14)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    
    # Add circle showing r0_turb = 0.8
    circle = plt.Circle((0, 0), 0.8, color='white', fill=False, linewidth=2, linestyle='--')
    ax1.add_patch(circle)
    ax1.text(0.85, 0.85, 'r₀=0.8', color='white', fontsize=10, 
             bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.7))
    
    plt.colorbar(im1, ax=ax1, label='|F|')
    
    # Panel 2: Force magnitude at x=0 (yz plane)
    ax2 = plt.subplot(3, 3, 2)
    i_slice = nx // 2
    fmag_yz = force_mag[:, :, i_slice]
    im2 = ax2.imshow(fmag_yz, origin='lower', cmap='viridis', 
                     extent=[-1, 1, -1, 1], interpolation='bilinear')
    ax2.set_title('Force Magnitude (x=0)', fontsize=14)
    ax2.set_xlabel('y')
    ax2.set_ylabel('z')
    
    circle2 = plt.Circle((0, 0), 0.8, color='white', fill=False, linewidth=2, linestyle='--')
    ax2.add_patch(circle2)
    
    plt.colorbar(im2, ax=ax2, label='|F|')
    
    # Panel 3: Force magnitude at y=0 (xz plane)
    ax3 = plt.subplot(3, 3, 3)
    j_slice = ny // 2
    fmag_xz = force_mag[:, j_slice, :]
    im3 = ax3.imshow(fmag_xz, origin='lower', cmap='viridis', 
                     extent=[-1, 1, -1, 1], interpolation='bilinear')
    ax3.set_title('Force Magnitude (y=0)', fontsize=14)
    ax3.set_xlabel('x')
    ax3.set_ylabel('z')
    
    circle3 = plt.Circle((0, 0), 0.8, color='white', fill=False, linewidth=2, linestyle='--')
    ax3.add_patch(circle3)
    
    plt.colorbar(im3, ax=ax3, label='|F|')
    
    # Panel 4: Radial profile along x-axis
    ax4 = plt.subplot(3, 3, 4)
    j_center = ny // 2
    k_center = nz // 2
    x_coords = x
    
    # Extract force components along x-axis
    fx_line = force_x[k_center, j_center, :]
    fy_line = force_y[k_center, j_center, :]
    fz_line = force_z[k_center, j_center, :]
    fmag_line = force_mag[k_center, j_center, :]
    
    ax4.plot(x_coords, fx_line, 'b-', linewidth=2, label='Fx')
    ax4.plot(x_coords, fy_line, 'r-', linewidth=2, label='Fy')
    ax4.plot(x_coords, fz_line, 'g-', linewidth=2, label='Fz')
    ax4.plot(x_coords, fmag_line, 'k--', linewidth=1.5, label='|F|')
    
    # Mark r0_turb boundaries
    ax4.axvline(x=-0.8, color='gray', linestyle='--', alpha=0.7)
    ax4.axvline(x=0.8, color='gray', linestyle='--', alpha=0.7)
    ax4.text(0.82, ax4.get_ylim()[1]*0.9, 'r₀', fontsize=10, color='gray')
    
    ax4.set_title('Force Components Along x-axis (y=z=0)', fontsize=14)
    ax4.set_xlabel('x')
    ax4.set_ylabel('Force')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    # Panel 5: Radial average of force magnitude
    ax5 = plt.subplot(3, 3, 5)
    
    # Compute radial distance for each point
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    R = np.sqrt(X**2 + Y**2 + Z**2)
    
    # Bin the data by radius
    r_bins = np.linspace(0, 1.5, 50)
    r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])
    force_avg = np.zeros_like(r_centers)
    
    for i in range(len(r_centers)):
        mask = (R >= r_bins[i]) & (R < r_bins[i+1])
        if np.sum(mask) > 0:
            force_avg[i] = np.mean(force_mag[mask.T])  # Transpose for correct indexing
    
    ax5.plot(r_centers, force_avg, 'k-', linewidth=2)
    ax5.axvline(x=0.8, color='red', linestyle='--', alpha=0.7, label='r₀=0.8')
    ax5.set_title('Radially Averaged Force Magnitude', fontsize=14)
    ax5.set_xlabel('r')
    ax5.set_ylabel('⟨|F|⟩')
    ax5.grid(True, alpha=0.3)
    ax5.legend()
    
    # Panel 6: Angular structure - spherical harmonic decomposition at fixed r
    ax6 = plt.subplot(3, 3, 6)
    
    # Extract shell at r ≈ 0.5
    r_target = 0.5
    shell_mask = (R >= r_target - 0.05) & (R <= r_target + 0.05)
    
    # Get theta and phi for points in shell
    theta = np.arccos(np.clip(Z/np.maximum(R, 1e-10), -1, 1))
    phi = np.arctan2(Y, X)
    
    # Simple angular power spectrum (l-modes)
    l_max = 10
    angular_power = []
    
    for l in range(l_max + 1):
        # Very simplified: just look at variation with theta
        theta_bins = np.linspace(0, np.pi, 20)
        theta_variation = []
        
        for i in range(len(theta_bins)-1):
            theta_mask = (theta >= theta_bins[i]) & (theta < theta_bins[i+1]) & shell_mask.T
            if np.sum(theta_mask) > 0:
                theta_variation.append(np.std(force_mag[theta_mask]))
        
        if theta_variation:
            angular_power.append(np.mean(theta_variation))
        else:
            angular_power.append(0)
    
    ax6.bar(range(l_max + 1), angular_power, color='steelblue', alpha=0.7)
    ax6.set_title(f'Angular Structure at r≈{r_target}', fontsize=14)
    ax6.set_xlabel('Angular scale (l-mode proxy)')
    ax6.set_ylabel('Force variation')
    ax6.grid(True, alpha=0.3, axis='y')
    
    # Panel 7-9: Vector field visualizations
    # Panel 7: Vector field in xy plane
    ax7 = plt.subplot(3, 3, 7)
    skip = max(4, nx//20)  # Adaptive skip for clarity
    
    X_2d, Y_2d = np.meshgrid(x[::skip], y[::skip])
    fx_slice = force_x[k_slice, ::skip, ::skip]
    fy_slice = force_y[k_slice, ::skip, ::skip]
    
    # Background with magnitude
    im7 = ax7.imshow(fmag_xy, origin='lower', cmap='gray', alpha=0.5,
                     extent=[-1, 1, -1, 1], interpolation='bilinear')
    
    # Quiver plot
    Q = ax7.quiver(X_2d, Y_2d, fx_slice, fy_slice,
                   force_mag[k_slice, ::skip, ::skip], cmap='plasma', 
                   scale=None, width=0.003)
    
    ax7.set_title('Force Vectors in xy-plane (z=0)', fontsize=14)
    ax7.set_xlabel('x')
    ax7.set_ylabel('y')
    ax7.set_aspect('equal')
    
    circle7 = plt.Circle((0, 0), 0.8, color='white', fill=False, linewidth=1, linestyle='--')
    ax7.add_patch(circle7)
    
    # Panel 8: Force components in spherical coordinates (r-component)
    ax8 = plt.subplot(3, 3, 8)
    
    # Compute radial component Fr = F·r̂
    R_safe = np.maximum(R, 1e-10)
    Fr = (force_x.T * X/R_safe + force_y.T * Y/R_safe + force_z.T * Z/R_safe).T
    
    Fr_xy = Fr[k_slice, :, :]
    im8 = ax8.imshow(Fr_xy, origin='lower', cmap='RdBu_r', 
                     extent=[-1, 1, -1, 1], interpolation='bilinear',
                     vmin=-np.max(np.abs(Fr_xy)), vmax=np.max(np.abs(Fr_xy)))
    ax8.set_title('Radial Force Component (z=0)', fontsize=14)
    ax8.set_xlabel('x')
    ax8.set_ylabel('y')
    
    circle8 = plt.Circle((0, 0), 0.8, color='black', fill=False, linewidth=2, linestyle='--')
    ax8.add_patch(circle8)
    
    plt.colorbar(im8, ax=ax8, label='Fr')
    
    # Panel 9: 2D FFT to show k-space structure
    ax9 = plt.subplot(3, 3, 9)
    
    # Compute 2D FFT of the force magnitude in xy plane
    fft_data = np.fft.fft2(fmag_xy)
    fft_data = np.fft.fftshift(fft_data)
    power_spectrum = np.abs(fft_data)**2
    
    # Create k-space coordinates
    kx = np.fft.fftfreq(nx, d=2.0/nx)
    ky = np.fft.fftfreq(ny, d=2.0/ny)
    kx = np.fft.fftshift(kx) * nx/2
    ky = np.fft.fftshift(ky) * ny/2
    
    # Plot power spectrum (log scale)
    im9 = ax9.imshow(np.log10(power_spectrum + 1e-10), origin='lower', cmap='hot',
                     extent=[kx.min(), kx.max(), ky.min(), ky.max()],
                     interpolation='bilinear', vmin=-5, vmax=5)
    ax9.set_title('Power Spectrum (log scale)', fontsize=14)
    ax9.set_xlabel('kx')
    ax9.set_ylabel('ky')
    ax9.set_xlim(-20, 20)
    ax9.set_ylim(-20, 20)
    
    plt.colorbar(im9, ax=ax9, label='log10(Power)')
    
    # Overall title with metadata
    time = data.get('Time', 0.0)
    cycle = data.get('NumCycles', 0)
    fig.suptitle(f'SFB Turbulence Driving Field Analysis\nTime = {time:.3f}, Cycle = {cycle}', 
                 fontsize=16)
    
    plt.tight_layout()
    plt.savefig(output_name, dpi=150, bbox_inches='tight')
    print(f"Saved visualization to {output_name}")
    
    # Print statistics
    print(f"\nStatistics:")
    print(f"  Max force magnitude: {force_mag.max():.3e}")
    print(f"  Mean force magnitude: {force_mag.mean():.3e}")
    print(f"  Force inside r₀: {force_mag[R.T < 0.8].mean():.3e}")
    print(f"  Force outside r₀: {force_mag[R.T >= 0.8].mean():.3e}")

def main():
    if len(sys.argv) > 1:
        bin_file = sys.argv[1]
    else:
        # Default to the latest output
        bin_files = sorted([f for f in os.listdir('.') if 'turb_force' in f and f.endswith('.bin')])
        if bin_files:
            bin_file = bin_files[-1]
        else:
            print("Error: No turbulence force binary files found!")
            return
    
    if not os.path.exists(bin_file):
        print(f"Error: File {bin_file} not found!")
        return
    
    output_name = bin_file.replace('.bin', '_sfb_viz.png')
    visualize_sfb_turb_field(bin_file, output_name)

if __name__ == "__main__":
    main()