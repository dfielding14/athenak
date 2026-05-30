# A simple script for converting collections of .bin/.cbin files to
# .athdf/.xdmf files using bin_convert.

# Python modules
import argparse
import os
import sys

# AthenaK modules
if __package__:
    from . import bin_convert
    from .io_reader_common import (
        ReaderLimits,
        normalize_reader_limits,
        require_live_bytes,
    )
else:
    import bin_convert
    from io_reader_common import ReaderLimits, normalize_reader_limits, require_live_bytes


def _bounded_batch_files(file_stem, limits=None):
    """Return deterministic .bin/.cbin candidates under the reader live budget."""
    limits = normalize_reader_limits(limits)
    directory, basename = os.path.split(file_stem)
    scan_directory = directory or "."
    files = []
    retained_member_bytes = 0
    retained_bytes = sys.getsizeof(files) + sys.getsizeof(scan_directory)
    require_live_bytes(retained_bytes, "batch binary inventory", limits)
    with os.scandir(scan_directory) as entries:
        for entry in entries:
            if (
                not entry.is_file()
                or not entry.name.startswith(basename)
                or not entry.name.endswith((".bin", ".cbin"))
                or entry.name.endswith(".sph.bin")
            ):
                continue
            path = os.path.join(directory, entry.name) if directory else entry.name
            item_bytes = sys.getsizeof(path)
            require_live_bytes(
                retained_bytes + item_bytes + 2 * sys.getsizeof(None),
                "batch binary inventory",
                limits,
            )
            files.append(path)
            retained_member_bytes += item_bytes
            retained_bytes = (
                retained_member_bytes
                + sys.getsizeof(files)
                + sys.getsizeof(scan_directory)
            )
            require_live_bytes(retained_bytes, "batch binary inventory", limits)
    require_live_bytes(
        retained_bytes + len(files) * sys.getsizeof(None),
        "batch binary inventory sort peak",
        limits,
    )
    files.sort()
    return files, retained_bytes


# Main function
def main(**kwargs):
    # Get the root name for the file.
    limits = normalize_reader_limits(kwargs.get('limits'))
    files, retained_bytes = _bounded_batch_files(kwargs['file_stem'], limits)
    if len(files) < 1:
        print(f"No files found with stem {kwargs['file_stem']}")
        quit()
    require_live_bytes(
        retained_bytes + 1, "batch binary conversion budget", limits
    )
    conversion_limits = ReaderLimits(
        max_live_bytes=limits.max_live_bytes - retained_bytes,
        max_header_read_bytes=limits.max_header_read_bytes,
        max_payload_read_bytes=limits.max_payload_read_bytes,
    )

    total = len(files)
    count = 1

    for fname in files:
        bin_convert.convert_file(
            fname,
            assemble_shards=kwargs.get('assemble_shards', False),
            limits=conversion_limits,
        )
        if kwargs['verbose']:
            print(f'Converting {count}/{total}: {fname}')
        count = count + 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'file_stem', help='path to files, excluding sequence number and .bin/.cbin suffix'
    )
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print file conversion progress')
    parser.add_argument(
        '--assemble-shards',
        action='store_true',
        help='assemble matching rank_* or node_* sibling shards before conversion',
    )
    bin_convert.add_reader_limit_arguments(parser)
    args = parser.parse_args()
    options = vars(args)
    options['limits'] = bin_convert.reader_limits_from_args(args)
    main(**options)
