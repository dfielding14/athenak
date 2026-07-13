#!/usr/bin/env python3
"""Hash a written GOTHAM conversion and atomically mark it complete."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import os
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return path, digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description='Hash and finalize one written GOTHAM restart conversion.')
    parser.add_argument('directory', type=Path, help='conversion output directory')
    parser.add_argument('--workers', type=int, default=16,
                        help='parallel hash workers (default: 16)')
    args = parser.parse_args()

    directory = args.directory.resolve()
    manifest_path = directory / 'conversion_manifest.txt'
    if not (directory / 'CONVERSION_WRITTEN').is_file() or not manifest_path.is_file():
        raise SystemExit('conversion has not reached its written state')

    manifest_text = manifest_path.read_text()
    metadata = dict(line.split('=', 1) for line in manifest_text.splitlines()
                    if '=' in line)
    expected_shards = int(metadata['output_shards'])
    fluid = list((directory / 'rst').glob('*.rst'))
    particles = list((directory / 'rst_prtcl').glob('*.rst_prtcl'))
    if len(fluid) != 1 or len(particles) != 1:
        raise SystemExit('expected exactly one fluid manifest and one particle restart')
    shards = sorted((directory / 'rst').glob('node_*/' + fluid[0].name))
    if len(shards) != expected_shards:
        raise SystemExit(f'expected {expected_shards} fluid shards, found {len(shards)}')

    files = fluid + particles + shards
    if any(path.stat().st_size == 0 for path in files):
        raise SystemExit('conversion contains an empty restart file')
    workers = max(1, min(args.workers, len(files)))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        hashes = list(pool.map(sha256, files))

    lines = [line for line in manifest_text.splitlines()
             if not line.startswith('hashes=') and not line.startswith('sha256.')]
    lines.append('hashes=sha256')
    for path, digest in hashes:
        lines.append(f'sha256.{path.relative_to(directory)}={digest}')

    manifest_tmp = directory / 'conversion_manifest.tmp'
    manifest_tmp.write_text('\n'.join(lines) + '\n')
    os.replace(manifest_tmp, manifest_path)
    marker_tmp = directory / 'CONVERSION_COMPLETE.tmp'
    marker_tmp.write_text('Restart pair, manifest, and SHA-256 inventory verified.\n')
    os.replace(marker_tmp, directory / 'CONVERSION_COMPLETE')


if __name__ == '__main__':
    main()
