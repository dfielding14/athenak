#!/usr/bin/env python3
"""
Visualize turbulence driving field as a 2D slice with refinement boundaries
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from bin_convert_new import read_binary_as_athdf

def visualize_turb_field(filename, output_name):
    """Create a nice visualization of the turbulence driving field"""
    
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
    
    # Take slice at z=nz//2
    k_slice = nz // 2
    fx_slice = force_x[k_slice, :, :]
    fy_slice = force_y[k_slice, :, :]
    fmag_slice = force_mag[k_slice, :, :]
    
    # Create figure with multiple panels
    fig = plt.figure(figsize=(16, 12))
    
    # Panel 1: Force magnitude with refinement regions
    ax1 = plt.subplot(2, 3, 1)
    im1 = ax1.imshow(fmag_slice, origin='lower', cmap='viridis', 
                     extent=[-0.5, 0.5, -0.5, 0.5], interpolation='bilinear')
    ax1.set_title('Turbulence Force Magnitude', fontsize=14)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    
    # Add refinement region boxes
    refinement_boxes = [
        patches.Rectangle((-0.3, -0.5), 0.2, 1.0, linewidth=2, 
                         edgecolor='white', facecolor='none', linestyle='--'),
        patches.Rectangle((0.1, -0.5), 0.2, 1.0, linewidth=2, 
                         edgecolor='white', facecolor='none', linestyle='--')
    ]
    for box in refinement_boxes:
        ax1.add_patch(box)
    ax1.text(-0.2, 0.45, 'Refined', color='white', ha='center', fontsize=10, 
             bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5))
    ax1.text(0.2, 0.45, 'Refined', color='white', ha='center', fontsize=10,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.5))
    
    plt.colorbar(im1, ax=ax1, label='|F|')
    
    # Panel 2: X-component of force
    ax2 = plt.subplot(2, 3, 2)
    im2 = ax2.imshow(fx_slice, origin='lower', cmap='RdBu_r', 
                     extent=[-0.5, 0.5, -0.5, 0.5], interpolation='bilinear',
                     vmin=-np.max(np.abs(fx_slice)), vmax=np.max(np.abs(fx_slice)))
    ax2.set_title('Force X-component', fontsize=14)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    for box in refinement_boxes:
        ax2.add_patch(patches.Rectangle(box.get_xy(), box.get_width(), box.get_height(),
                                       linewidth=2, edgecolor='black', facecolor='none', linestyle='--'))
    plt.colorbar(im2, ax=ax2, label='Fx')
    
    # Panel 3: Y-component of force
    ax3 = plt.subplot(2, 3, 3)
    im3 = ax3.imshow(fy_slice, origin='lower', cmap='RdBu_r', 
                     extent=[-0.5, 0.5, -0.5, 0.5], interpolation='bilinear',
                     vmin=-np.max(np.abs(fy_slice)), vmax=np.max(np.abs(fy_slice)))
    ax3.set_title('Force Y-component', fontsize=14)
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    for box in refinement_boxes:
        ax3.add_patch(patches.Rectangle(box.get_xy(), box.get_width(), box.get_height(),
                                       linewidth=2, edgecolor='black', facecolor='none', linestyle='--'))
    plt.colorbar(im3, ax=ax3, label='Fy')
    
    # Panel 4: Vector field (quiver plot)
    ax4 = plt.subplot(2, 3, 4)
    # Downsample for clarity
    skip = 4
    x = np.linspace(-0.5, 0.5, nx)
    y = np.linspace(-0.5, 0.5, ny)
    X, Y = np.meshgrid(x[::skip], y[::skip])
    
    # Background with magnitude
    im4 = ax4.imshow(fmag_slice, origin='lower', cmap='gray', alpha=0.5,
                     extent=[-0.5, 0.5, -0.5, 0.5], interpolation='bilinear')
    
    # Quiver plot
    Q = ax4.quiver(X, Y, fx_slice[::skip, ::skip], fy_slice[::skip, ::skip],
                   fmag_slice[::skip, ::skip], cmap='plasma', scale=2.0, width=0.003)
    ax4.set_title('Force Vector Field', fontsize=14)
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.set_aspect('equal')
    for box in refinement_boxes:
        ax4.add_patch(patches.Rectangle(box.get_xy(), box.get_width(), box.get_height(),
                                       linewidth=2, edgecolor='white', facecolor='none', linestyle='--'))
    plt.colorbar(Q, ax=ax4, label='|F|')
    
    # Panel 5: Line cut through y=0
    ax5 = plt.subplot(2, 3, 5)
    j_center = ny // 2
    x_coords = np.linspace(-0.5, 0.5, nx)
    
    ax5.plot(x_coords, fx_slice[j_center, :], 'b-', linewidth=2, label='Fx')
    ax5.plot(x_coords, fy_slice[j_center, :], 'r-', linewidth=2, label='Fy')
    ax5.plot(x_coords, fmag_slice[j_center, :], 'k--', linewidth=1.5, label='|F|')
    
    ax5.set_title('Force Components at y=0', fontsize=14)
    ax5.set_xlabel('x')
    ax5.set_ylabel('Force')
    ax5.grid(True, alpha=0.3)
    ax5.legend()
    
    # Mark refinement boundaries
    for x_bound in [-0.3, -0.1, 0.1, 0.3]:
        ax5.axvline(x=x_bound, color='gray', linestyle='--', alpha=0.5)
    
    # Panel 6: Fourier spectrum
    ax6 = plt.subplot(2, 3, 6)
    
    # Compute 2D FFT of the force magnitude
    fft_data = np.fft.fft2(fmag_slice)
    fft_data = np.fft.fftshift(fft_data)
    power_spectrum = np.abs(fft_data)**2
    
    # Create k-space coordinates
    kx = np.fft.fftfreq(nx, d=1.0/nx)
    ky = np.fft.fftfreq(ny, d=1.0/ny)
    kx = np.fft.fftshift(kx) * nx
    ky = np.fft.fftshift(ky) * ny
    
    # Plot power spectrum (log scale)
    im6 = ax6.imshow(np.log10(power_spectrum + 1e-10), origin='lower', cmap='hot',
                     extent=[kx.min(), kx.max(), ky.min(), ky.max()],
                     interpolation='bilinear', vmin=-5, vmax=5)
    ax6.set_title('Power Spectrum (log scale)', fontsize=14)
    ax6.set_xlabel('kx')
    ax6.set_ylabel('ky')
    ax6.set_xlim(-10, 10)
    ax6.set_ylim(-10, 10)
    
    # Mark the driven modes
    circle_driven = plt.Circle((0, 0), 4, color='cyan', fill=False, linewidth=2, linestyle='--')
    ax6.add_patch(circle_driven)
    ax6.text(0, -8, 'Driven modes\n(2 ≤ |k| ≤ 4)', ha='center', color='cyan',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.7))
    
    plt.colorbar(im6, ax=ax6, label='log10(Power)')
    
    # Overall title with metadata
    time = data.get('Time', 0.0)
    cycle = data.get('NumCycles', 0)
    fig.suptitle(f'Turbulence Driving Field Analysis\nTime = {time:.3f}, Cycle = {cycle}', 
                 fontsize=16)
    
    plt.tight_layout()
    plt.savefig(output_name, dpi=150, bbox_inches='tight')
    print(f"Saved visualization to {output_name}")
    
    # Print some statistics
    print(f"\nStatistics:")
    print(f"  Max force magnitude: {fmag_slice.max():.3e}")
    print(f"  Mean force magnitude: {fmag_slice.mean():.3e}")
    print(f"  Force components range:")
    print(f"    Fx: [{fx_slice.min():.3e}, {fx_slice.max():.3e}]")
    print(f"    Fy: [{fy_slice.min():.3e}, {fy_slice.max():.3e}]")

def main():
    if len(sys.argv) > 1:
        bin_file = sys.argv[1]
    else:
        # Default to the latest output
        bin_files = sorted([f for f in os.listdir('bin') if 'turb_force' in f and f.endswith('.bin')])
        if bin_files:
            bin_file = os.path.join('bin', bin_files[-1])
        else:
            print("Error: No turbulence force binary files found!")
            return
    
    if not os.path.exists(bin_file):
        print(f"Error: File {bin_file} not found!")
        return
    
    output_name = bin_file.replace('.bin', '_field_viz.png').replace('bin/', '')
    visualize_turb_field(bin_file, output_name)

if __name__ == "__main__":
    main()