#!/usr/bin/env python3
"""
Analyze particle movement across VTK particle files.

This script reads particle VTK files (same reader as the plot script) and prints
per-step verbose positions for selected particle tags and a summary of net
displacement from fstart to fend. It was extracted from plot_particle_plane.py
so users can run analysis separately from plotting.
"""

import numpy as np
import argparse
import sys
from pathlib import Path

# Reuse the reader from the plotting script by importing it if available. If not,
# copy a minimal reader here to avoid circular imports.
try:
    from plot_particle_plane import read_vtk_particle_file
except Exception:
    # Minimal local copy of reader for standalone operation
    def read_vtk_particle_file(filename):
        with open(filename, 'rb') as f:
            header_lines = []
            while True:
                line = f.readline().decode('ascii')
                header_lines.append(line.strip())
                if 'POINTS' in line:
                    parts = line.split()
                    nparticles = int(parts[1])
                    break
                if len(header_lines) > 100:
                    raise ValueError(f"Could not find POINTS line in {filename}")
            time = None
            cycle = None
            for line in header_lines:
                if 'time=' in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'time=' and i+1 < len(parts):
                            time = float(parts[i+1])
                        if part.startswith('cycle='):
                            cycle = int(part.split('=')[1])
            position_data = f.read(4 * 3 * nparticles)
            positions = np.frombuffer(position_data, dtype='>f4').reshape((nparticles, 3))
            # Skip to POINT_DATA
            while True:
                line = f.readline().decode('ascii')
                if 'POINT_DATA' in line:
                    break
                if not line:
                    raise ValueError(f"Could not find POINT_DATA in {filename}")
            gid = None
            ptag = None
            while True:
                line = f.readline()
                if not line:
                    break
                line = line.decode('ascii').strip()
                if 'SCALARS gid' in line:
                    f.readline()
                    gid_data = f.read(4 * nparticles)
                    gid = np.frombuffer(gid_data, dtype='>f4').astype(np.int32)
                elif 'SCALARS ptag' in line:
                    f.readline()
                    ptag_data = f.read(4 * nparticles)
                    ptag = np.frombuffer(ptag_data, dtype='>f4').astype(np.int32)
        return {'time': time, 'cycle': cycle, 'nparticles': nparticles, 'positions': positions, 'gid': gid, 'ptag': ptag}


def analyze(directories, fstart, fend, basename, ptag_filter, verbose, verbose_summary, displacement_threshold, nparticles):
    if isinstance(directories, (str, Path)):
        directories = [directories]

    for directory in directories:
        directory = Path(directory)
        if not directory.exists():
            print(f"Warning: Directory {directory} does not exist, skipping")
            continue

        print(f"\nAnalyzing directory: {directory}")

        first_positions = {}
        last_positions = {}
        
        # If nparticles is specified but no ptag_filter, randomly sample from first file
        random_sample_tags = None
        if nparticles is not None and ptag_filter is None:
            first_file = directory / f"{basename}.{fstart:05d}.part.vtk"
            if first_file.exists():
                first_data = read_vtk_particle_file(first_file)
                if first_data['ptag'] is not None:
                    available_tags = first_data['ptag']
                    n_available = len(available_tags)
                    n_sample = min(nparticles, n_available)
                    # Randomly sample without replacement
                    np.random.seed(42)  # For reproducibility
                    sample_indices = np.random.choice(n_available, size=n_sample, replace=False)
                    random_sample_tags = available_tags[sample_indices]
                    ptag_filter = random_sample_tags.tolist()
                    print(f"  Randomly sampled {n_sample} particles from {n_available} total")
                    print(f"  Selected PTAGs: {sorted(ptag_filter)[:10]}{'...' if len(ptag_filter) > 10 else ''}")
                else:
                    print(f"  Warning: No PTAG data found in {first_file}")

        for step_num in range(fstart, fend + 1):
            filename = directory / f"{basename}.{step_num:05d}.part.vtk"
            if not filename.exists():
                print(f"  Warning: {filename} not found, stopping at step {step_num-1}")
                break

            data = read_vtk_particle_file(filename)
            if ptag_filter is not None and data['ptag'] is not None:
                mask = np.isin(data['ptag'], ptag_filter)
                positions = data['positions'][mask]
                filtered_ptags = data['ptag'][mask]
                n_filtered = np.sum(mask)
                print(f"  Step {step_num}: time={data['time']:.6e}, {n_filtered}/{data['nparticles']} particles (filtered)")

                if verbose and n_filtered > 0:
                    print(f"    Tagged particle positions:")
                    sort_idx = np.argsort(filtered_ptags)
                    for idx in sort_idx:
                        ptag = filtered_ptags[idx]
                        pos = positions[idx]
                        print(f"      PTAG {ptag:6d}: x={pos[0]:18.15f}, y={pos[1]:18.15f}, z={pos[2]:18.15f}")

                if n_filtered > 0:
                    for idx, ptag in enumerate(filtered_ptags):
                        if ptag not in first_positions:
                            first_positions[ptag] = (positions[idx].copy(), data['time'])
                        last_positions[ptag] = (positions[idx].copy(), data['time'])
            else:
                print(f"  Step {step_num}: time={data['time']:.6e}, {data['nparticles']} particles")

        # Summary
        if verbose_summary and len(first_positions) > 0:
            print(f"\n=== Particle Displacement Summary for {directory} ===")
            print(f"Displacement threshold: {displacement_threshold:.6f}")
            print(f"{'PTAG':>6s} | {'Time Start':>14s} | {'Time End':>14s} | {'Delta t':>14s} | {'dx':>18s} | {'dy':>18s} | {'dz':>18s} | {'|disp|':>18s} | {'Moved?':>7s}")
            particles_moved = 0
            total_particles = len(first_positions)
            for ptag in sorted(first_positions.keys()):
                pos_start, time_start = first_positions[ptag]
                pos_end, time_end = last_positions[ptag]
                dx = pos_end[0] - pos_start[0]
                dy = pos_end[1] - pos_start[1]
                dz = pos_end[2] - pos_start[2]
                displacement = np.sqrt(dx**2 + dy**2 + dz**2)
                dt = time_end - time_start
                moved = 'YES' if displacement >= displacement_threshold else 'NO'
                if displacement >= displacement_threshold:
                    particles_moved += 1
                print(f"{ptag:6d} | {time_start:14.6e} | {time_end:14.6e} | {dt:14.6e} | {dx:18.15f} | {dy:18.15f} | {dz:18.15f} | {displacement:18.15f} | {moved:>7s}")
            print()
            print(f"Total particles tracked: {total_particles}")
            print(f"Particles moved (|disp| >= {displacement_threshold:.6f}): {particles_moved} ({100*particles_moved/total_particles:.1f}%)")
            print(f"Particles stationary (|disp| < {displacement_threshold:.6f}): {total_particles - particles_moved} ({100*(total_particles - particles_moved)/total_particles:.1f}%)")


def main():
    parser = argparse.ArgumentParser(description='Analyze particle movement from pvtk files')
    parser.add_argument('--dir', type=str, nargs='+', required=True,
                        help='Directory(ies) containing VTK particle files')
    parser.add_argument('--fstart', type=int, required=True, help='First file number to read')
    parser.add_argument('--fend', type=int, required=True, help='Last file number to read')
    parser.add_argument('--basename', type=str, default='KH.prtcl_all', help='Base name of particle files')
    parser.add_argument('--ptag', type=int, nargs='+', help='Specific particle tag(s) to analyze')
    parser.add_argument('--nparticles', type=int, help='If --ptag not specified, randomly sample this many particles')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print per-step particle positions')
    parser.add_argument('--verbose-summary', action='store_true', help='Print summary of net displacement')
    parser.add_argument('--displacement-threshold', type=float, default=0.01, help='Threshold for counting particles as moved')

    args = parser.parse_args()

    if args.fstart > args.fend:
        print('Error: fstart must be <= fend')
        sys.exit(1)

    analyze(directories=args.dir, fstart=args.fstart, fend=args.fend, basename=args.basename,
            ptag_filter=args.ptag, verbose=args.verbose, verbose_summary=args.verbose_summary,
            displacement_threshold=args.displacement_threshold, nparticles=args.nparticles)


if __name__ == '__main__':
    sys.exit(main())
