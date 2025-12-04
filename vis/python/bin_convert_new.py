"""bin_convert_new.py
Utility helpers for working with Athena++ binary output files.

This module offers a high-performance pure-Python reader for Athena++ *.bin
snapshots, supporting:
  - Mesh-refined outputs (multiple AMR levels).
  - Coarsened (on-the-fly averaged) outputs.
  - Multi-rank dumps written with one file per MPI rank.
  - On-the-fly translation to an athdf-like in-memory representation.
  - Generation of companion XDMF files for ParaView/VisIt visualisation.

The original reader was developed by Lev Arzamasskiy (2021-11-15) and later
extended by George Wong (2022-01-27) and Drummond Fielding (2024-09-09).
This docstring and accompanying clean-up were added in 2025-06-19 to improve
clarity, bring the file closer to PEP-8 compliance, and document the public
API.  Down-stream code should rely only on the symbols exposed via ``__all__``
(defined at the bottom of the file).
"""
import numpy as np
import os
import glob
import io
import h5py
from typing import Any, Callable, Dict, List, Optional, Sequence

# --------------------------------------------------------------------------------------
# Module-level constants ----------------------------------------------------------------
# --------------------------------------------------------------------------------------

CODE_HEADER_EXPECTED: bytes = b"Athena"
SUPPORTED_VERSION: bytes = b"1.1"

# --------------------------------------------------------------------------------------
# Private utility helpers ---------------------------------------------------------------
# --------------------------------------------------------------------------------------


def _parse_header_to_dict(header: List[str]) -> Dict[str, Dict[str, str]]:
    """Parse header lines into a nested dictionary for O(1) lookups.

    Returns a dict mapping block names (e.g. "<mesh>") to dicts of key-value pairs.
    """
    data: Dict[str, Dict[str, str]] = {}
    current_block = "<none>"
    for line in header:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("<"):
            current_block = line
            continue
        if "=" in line:
            k, v = line.split("=", 1)
            if current_block not in data:
                data[current_block] = {}
            data[current_block][k.strip()] = v.strip()
    return data


def _get_from_header(header_dict: Dict[str, Dict[str, str]], blockname: str,
                     keyname: str) -> str:
    """Return the value for *keyname* inside *blockname* using O(1) lookup.

    Parameters
    ----------
    header_dict
        Parsed header dictionary from _parse_header_to_dict.
    blockname
        Name of the block (e.g. ``"<mesh>"``). Leading ``<`` and trailing ``>``
        are optional.
    keyname
        Parameter to extract (e.g. ``"nx1"``).
    """
    blockname = blockname.strip()
    if not blockname.startswith("<"):
        blockname = "<" + blockname
    if not blockname.endswith(">"):
        blockname += ">"

    try:
        return header_dict[blockname][keyname]
    except KeyError:
        raise KeyError(f"no parameter called {blockname}/{keyname}")


# --------------------------------------------------------------------------------------
# Shared low-level binary reader --------------------------------------------------------
# --------------------------------------------------------------------------------------


def _read_binary_impl(filename: str, *, coarsened: bool = False) -> Dict[str, Any]:
    """Parse a single Athena++ ``*.bin`` dump.

    This consolidates the duplicated logic previously found in
    ``read_binary`` and ``read_coarsened_binary``.  Both public wrappers now
    delegate to this function so that format tweaks touch just one place.
    """

    # Load file into an in-memory buffer for fast seeks
    with open(filename, "rb") as fh:
        raw_bytes = fh.read()

    fp = io.BytesIO(raw_bytes)
    filesize = len(raw_bytes)

    # ----------------------------- Header tokens -----------------------------
    tokens = fp.readline().split()
    if not tokens or tokens[0] != CODE_HEADER_EXPECTED:
        raise TypeError("Not an Athena++ binary dump (missing magic header)")

    version = tokens[-1].split(b"=")[-1]
    if version != SUPPORTED_VERSION:
        raise TypeError(
            f"Unsupported Athena++ binary version {version.decode()} "
            f"(expected {SUPPORTED_VERSION.decode()})"
        )

    # Parameter header lines
    pheader_lines = int(fp.readline().split(b"=")[-1])
    pheader: Dict[str, str] = {}
    for _ in range(pheader_lines - 1):
        k, v = [s.strip() for s in fp.readline().decode().split("=")]
        pheader[k] = v

    # Scalar meta-data
    time = float(pheader["time"])
    cycle = int(pheader["cycle"])
    locsizebytes = int(pheader["size of location"])
    varsizebytes = int(pheader["size of variable"])
    coarsen_factor = int(pheader.get("coarsening factor", "1"))

    if coarsened and coarsen_factor == 1:
        raise ValueError(
            "Reader asked for a coarsened dump but header lacks coarsening factor"
        )

    # Variable list
    nvars = int(fp.readline().split(b"=")[-1])
    var_names: List[str] = [tok.decode() for tok in fp.readline().split()[1:]]

    # Simulation header (ASCII)
    ascii_header_size = int(fp.readline().split(b"=")[-1])
    header_bytes = fp.read(ascii_header_size)
    header_lines = [
        ln.decode().split("#")[0].strip()
        for ln in header_bytes.split(b"\n") if ln
    ]

    # Parse header once for O(1) lookups
    header_dict = _parse_header_to_dict(header_lines)

    # Validate float sizes
    if locsizebytes not in (4, 8) or varsizebytes not in (4, 8):
        raise ValueError("Unsupported float byte sizes in header")

    dtype_loc = np.float64 if locsizebytes == 8 else np.float32
    dtype_var = np.float64 if varsizebytes == 8 else np.float32

    # Mesh sizes
    Nx1 = int(_get_from_header(header_dict, "<mesh>", "nx1"))
    Nx2 = int(_get_from_header(header_dict, "<mesh>", "nx2"))
    Nx3 = int(_get_from_header(header_dict, "<mesh>", "nx3"))
    nx1 = int(_get_from_header(header_dict, "<meshblock>", "nx1"))
    nx2 = int(_get_from_header(header_dict, "<meshblock>", "nx2"))
    nx3 = int(_get_from_header(header_dict, "<meshblock>", "nx3"))
    nghost = int(_get_from_header(header_dict, "<mesh>", "nghost"))

    x1min = float(_get_from_header(header_dict, "<mesh>", "x1min"))
    x1max = float(_get_from_header(header_dict, "<mesh>", "x1max"))
    x2min = float(_get_from_header(header_dict, "<mesh>", "x2min"))
    x2max = float(_get_from_header(header_dict, "<mesh>", "x2max"))
    x3min = float(_get_from_header(header_dict, "<mesh>", "x3min"))
    x3max = float(_get_from_header(header_dict, "<mesh>", "x3max"))

    # ----------------------------- Mesh blocks -----------------------------
    mb_index, mb_logical, mb_geometry = [], [], []
    mb_data: Dict[str, List[np.ndarray]] = {v: [] for v in var_names}

    while fp.tell() < filesize:
        mb_index.append(
            np.frombuffer(fp.read(24), dtype=np.int32).astype(np.int64) - nghost
        )
        nx1_out = (mb_index[-1][1] - mb_index[-1][0]) + 1
        nx2_out = (mb_index[-1][3] - mb_index[-1][2]) + 1
        nx3_out = (mb_index[-1][5] - mb_index[-1][4]) + 1

        mb_logical.append(np.frombuffer(fp.read(16), dtype=np.int32))
        mb_geometry.append(np.frombuffer(fp.read(6 * locsizebytes), dtype=dtype_loc))

        # Cell data block
        cells = nx1_out * nx2_out * nx3_out * nvars
        bytes_needed = cells * np.dtype(dtype_var).itemsize
        block_bytes = fp.read(bytes_needed)
        if len(block_bytes) != bytes_needed:
            raise IOError(
                "Unexpected EOF while reading mesh-block data; "
                "binary dump may be corrupted"
            )

        block = np.frombuffer(block_bytes, dtype=dtype_var, count=cells).reshape(
            nvars, nx3_out, nx2_out, nx1_out
        )
        for vi, v in enumerate(var_names):
            mb_data[v].append(block[vi])

    # Coarsening factor handling
    factor = coarsen_factor if coarsened else 1

    result: Dict[str, Any] = {
        "header": header_lines,
        "time": time,
        "cycle": cycle,
        "var_names": var_names,
        "nvars": nvars,
        "Nx1": Nx1 // factor,
        "Nx2": Nx2 // factor,
        "Nx3": Nx3 // factor,
        "x1min": x1min,
        "x1max": x1max,
        "x2min": x2min,
        "x2max": x2max,
        "x3min": x3min,
        "x3max": x3max,
        "n_mbs": len(mb_index),
        "nx1_mb": nx1 // factor,
        "nx2_mb": nx2 // factor,
        "nx3_mb": nx3 // factor,
        "nx1_out_mb": (mb_index[0][1] - mb_index[0][0]) + 1 if mb_index else 0,
        "nx2_out_mb": (mb_index[0][3] - mb_index[0][2]) + 1 if mb_index else 0,
        "nx3_out_mb": (mb_index[0][5] - mb_index[0][4]) + 1 if mb_index else 0,
        "mb_index": np.array(mb_index),
        "mb_logical": np.array(mb_logical),
        "mb_geometry": np.array(mb_geometry),
        "mb_data": mb_data,
    }
    if coarsened:
        result["number_of_moments"] = int(pheader["number of moments"])
    return result


# --------------------------------------------------------------------------------------
# Public binary readers -----------------------------------------------------------------
# --------------------------------------------------------------------------------------


def read_binary(filename: str) -> Dict[str, Any]:
    """Parse an Athena++ ``*.bin`` file produced **without** on-the-fly
    coarsening and return a dictionary identical to what `athdf` would
    contain on disk.
    """
    return _read_binary_impl(filename, coarsened=False)


def read_coarsened_binary(filename: str) -> Dict[str, Any]:
    """Parse a *coarsened* Athena++ ``*.bin`` snapshot produced with
    ``output_coarsening`` enabled."""
    return _read_binary_impl(filename, coarsened=True)


def _combine_rank_files(rank_files: List[str], reader_func: Callable) -> Dict[str, Any]:
    """Combine binary files from multiple ranks into a single dictionary.

    Parameters
    ----------
    rank_files
        Sorted list of rank file paths.
    reader_func
        Function to read a single binary file (read_binary or read_coarsened_binary).

    Returns
    -------
    combined_filedata
        Dictionary of combined fluid file data from all ranks.
    """
    # Read rank 0 to get metadata
    rank0_filedata = reader_func(rank_files[0])

    # Initialize with rank 0 data (avoid re-reading it)
    combined_filedata = rank0_filedata.copy()
    combined_filedata["mb_index"] = list(rank0_filedata["mb_index"])
    combined_filedata["mb_logical"] = list(rank0_filedata["mb_logical"])
    combined_filedata["mb_geometry"] = list(rank0_filedata["mb_geometry"])
    combined_filedata["mb_data"] = {
        var: list(rank0_filedata["mb_data"][var])
        for var in rank0_filedata["var_names"]
    }

    # Read remaining ranks (skip rank 0)
    for rank_filename in rank_files[1:]:
        rank_filedata = reader_func(rank_filename)
        combined_filedata["mb_index"].extend(rank_filedata["mb_index"])
        combined_filedata["mb_logical"].extend(rank_filedata["mb_logical"])
        combined_filedata["mb_geometry"].extend(rank_filedata["mb_geometry"])
        for var in rank0_filedata["var_names"]:
            combined_filedata["mb_data"][var].extend(rank_filedata["mb_data"][var])

    # Convert lists to numpy arrays
    combined_filedata["mb_index"] = np.array(combined_filedata["mb_index"])
    combined_filedata["mb_logical"] = np.array(combined_filedata["mb_logical"])
    combined_filedata["mb_geometry"] = np.array(combined_filedata["mb_geometry"])
    for var in rank0_filedata["var_names"]:
        combined_filedata["mb_data"][var] = np.array(
            combined_filedata["mb_data"][var]
        )

    combined_filedata["n_mbs"] = len(combined_filedata["mb_index"])
    return combined_filedata


def _find_rank_files(rank0_filename: str) -> List[str]:
    """Find all rank files matching the rank0 pattern."""
    pattern = (
        os.path.dirname(rank0_filename).replace("rank_00000000", "rank_*")
        + "/" + os.path.basename(rank0_filename)
    )
    rank_files = sorted(glob.glob(pattern))

    # Filter out empty files (e.g., slice outputs where rank has no intersecting meshblocks)
    # We check if files have any meshblock data by comparing to a minimum header size
    # A file with just headers and no meshblocks will be small (~few KB)
    MIN_DATA_SIZE = 1024  # Files smaller than this likely have no meshblock data
    file_sizes = np.array([os.path.getsize(f) for f in rank_files])
    small_files = file_sizes < MIN_DATA_SIZE
    if np.any(small_files) and not np.all(small_files):
        n_filtered = np.sum(small_files)
        rank_files = [f for f, is_small in zip(rank_files, small_files) if not is_small]
        print(f"Note: Filtered out {n_filtered} empty rank file(s) (slice output with no data)")

    return rank_files


def read_all_ranks_binary(rank0_filename: str) -> Dict[str, Any]:
    """Reads binary files from all ranks and combines them into a single dictionary.

    Parameters
    ----------
    rank0_filename
        Filename of the rank 0 binary file.

    Returns
    -------
    combined_filedata
        Dictionary of combined fluid file data from all ranks.
    """
    rank_files = _find_rank_files(rank0_filename)
    return _combine_rank_files(rank_files, read_binary)


def read_all_ranks_coarsened_binary(rank0_filename: str) -> Dict[str, Any]:
    """Reads coarsened binary files from all ranks and combines them.

    Parameters
    ----------
    rank0_filename
        Filename of the rank 0 binary file.

    Returns
    -------
    combined_filedata
        Dictionary of combined fluid file data from all ranks.
    """
    rank_files = _find_rank_files(rank0_filename)
    return _combine_rank_files(rank_files, read_coarsened_binary)


# --------------------------------------------------------------------------------------
# Shared athdf conversion logic ---------------------------------------------------------
# --------------------------------------------------------------------------------------


def _organize_as_athdf(
    filedata: Dict[str, Any],
    *,
    data: Optional[Dict[str, Any]] = None,
    quantities: Optional[Sequence[str]] = None,
    dtype: Optional[np.dtype] = None,
    level: Optional[int] = None,
    return_levels: bool = False,
    subsample: bool = False,
    fast_restrict: bool = False,
    x1_min: Optional[float] = None,
    x1_max: Optional[float] = None,
    x2_min: Optional[float] = None,
    x2_max: Optional[float] = None,
    x3_min: Optional[float] = None,
    x3_max: Optional[float] = None,
    vol_func: Optional[Callable] = None,
    center_func_1: Optional[Callable] = None,
    center_func_2: Optional[Callable] = None,
    center_func_3: Optional[Callable] = None,
    num_ghost: int = 0,
    use_output_block_size: bool = True,
) -> Dict[str, Any]:
    """Convert raw binary filedata to athdf-like dictionary format.

    This is the shared implementation for all *_as_athdf functions.

    Parameters
    ----------
    filedata
        Raw binary data from read_binary or similar.
    use_output_block_size
        If True, use nx*_out_mb for block sizes (standard binary).
        If False, use nx*_mb (coarsened binary).
    """
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    if dtype is None:
        dtype = np.float32

    # Extract size information
    max_level = int(np.max(filedata['mb_logical'][:, 3]))
    if level is None:
        level = max_level

    if use_output_block_size:
        block_size = [
            filedata['nx1_out_mb'],
            filedata['nx2_out_mb'],
            filedata['nx3_out_mb']
        ]
    else:
        block_size = [filedata['nx1_mb'], filedata['nx2_mb'], filedata['nx3_mb']]

    root_grid_size = [filedata['Nx1'], filedata['Nx2'], filedata['Nx3']]
    levels = filedata['mb_logical'][:, 3]
    logical_locations = filedata['mb_logical'][:, :3]

    # Calculate output grid dimensions
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:  # sum or slice
            other_locations = list(zip(
                levels,
                logical_locations[:, (d + 1) % 3],
                logical_locations[:, (d + 2) % 3]
            ))
            if len(set(other_locations)) == len(other_locations):  # effective slice
                nx_vals.append(1)
            else:  # nontrivial sum
                num_blocks_this_dim = 0
                for level_this_dim, loc_this_dim in zip(levels, logical_locations[:, d]):
                    if level_this_dim <= level:
                        possible_max = (loc_this_dim + 1) * 2 ** (level - level_this_dim)
                    else:
                        possible_max = (loc_this_dim + 1) // 2 ** (level_this_dim - level)
                    num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                nx_vals.append(num_blocks_this_dim)
        elif block_size[d] == 1:  # singleton dimension
            nx_vals.append(1)
        else:  # normal case
            nx_vals.append(root_grid_size[d] * 2 ** level + 2 * num_ghost)

    nx1, nx2, nx3 = nx_vals
    lx1 = nx1 // block_size[0] if block_size[0] > 0 else 0
    lx2 = nx2 // block_size[1] if block_size[1] > 0 else 0
    lx3 = nx3 // block_size[2] if block_size[2] > 0 else 0

    # Default coordinate functions
    if vol_func is None:
        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)
    if center_func_1 is None:
        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)
    if center_func_2 is None:
        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)
    if center_func_3 is None:
        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    center_funcs = [center_func_1, center_func_2, center_func_3]

    # Populate coordinate arrays
    for d in range(1, 4):
        xf = f'x{d}f'
        xv = f'x{d}v'
        nx = nx_vals[d - 1]
        xmin = filedata[f'x{d}min']
        xmax = filedata[f'x{d}max']
        if nx == 1:
            data[xf] = np.array([xmin, xmax], dtype=dtype)
        else:
            data[xf] = np.linspace(xmin, xmax, nx + 1, dtype=dtype)
        data[xv] = np.empty(nx, dtype=dtype)
        for i in range(nx):
            data[xv][i] = center_funcs[d - 1](data[xf][i], data[xf][i + 1])

    # Create list of quantities
    if quantities is None:
        quantities = filedata['var_names']

    # Account for spatial selection
    i_min, i_max = 0, nx1
    j_min, j_max = 0, nx2
    k_min, k_max = 0, nx3
    if x1_min is not None:
        i_min = max(i_min, np.searchsorted(data['x1f'], x1_min))
    if x1_max is not None:
        i_max = min(i_max, np.searchsorted(data['x1f'], x1_max))
    if x2_min is not None:
        j_min = max(j_min, np.searchsorted(data['x2f'], x2_min))
    if x2_max is not None:
        j_max = min(j_max, np.searchsorted(data['x2f'], x2_max))
    if x3_min is not None:
        k_min = max(k_min, np.searchsorted(data['x3f'], x3_min))
    if x3_max is not None:
        k_max = min(k_max, np.searchsorted(data['x3f'], x3_max))

    # Prepare arrays for data and bookkeeping
    if new_data:
        for q in quantities:
            data[q] = np.zeros((k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype)
        if return_levels:
            data['Levels'] = np.empty(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)

    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Process each MeshBlock
    for block_num in range(filedata['n_mbs']):
        block_level = levels[block_num]
        block_location = logical_locations[block_num]

        if block_level <= level:
            # Prolongation case: block at same level or coarser
            s = 2 ** (level - block_level)
            il_d = block_location[0] * block_size[0] * s if nx1 > 1 else 0
            jl_d = block_location[1] * block_size[1] * s if nx2 > 1 else 0
            kl_d = block_location[2] * block_size[2] * s if nx3 > 1 else 0
            iu_d = il_d + block_size[0] * s if nx1 > 1 else 1
            ju_d = jl_d + block_size[1] * s if nx2 > 1 else 1
            ku_d = kl_d + block_size[2] * s if nx3 > 1 else 1

            il_s = max(il_d, i_min) - il_d
            jl_s = max(jl_d, j_min) - jl_d
            kl_s = max(kl_d, k_min) - kl_d
            iu_s = min(iu_d, i_max) - il_d
            ju_s = min(ju_d, j_max) - jl_d
            ku_s = min(ku_d, k_max) - kl_d

            if il_s >= iu_s or jl_s >= ju_s or kl_s >= ku_s:
                continue

            il_d = max(il_d, i_min) - i_min
            jl_d = max(jl_d, j_min) - j_min
            kl_d = max(kl_d, k_min) - k_min
            iu_d = min(iu_d, i_max) - i_min
            ju_d = min(ju_d, j_max) - j_min
            ku_d = min(ku_d, k_max) - k_min

            for q in quantities:
                block_data = filedata['mb_data'][q][block_num]
                if s > 1:
                    block_data = np.repeat(
                        np.repeat(
                            np.repeat(block_data, s, axis=2),
                            s, axis=1
                        ),
                        s, axis=0
                    )
                data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = \
                    block_data[kl_s:ku_s, jl_s:ju_s, il_s:iu_s]

            if return_levels:
                data['Levels'][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level

        else:
            # Restriction case: block is finer than requested level
            if subsample:
                # Simple subsampling: take every s-th cell
                s = 2 ** (block_level - level)
                il_d = block_location[0] // s * block_size[0] if nx1 > 1 else 0
                jl_d = block_location[1] // s * block_size[1] if nx2 > 1 else 0
                kl_d = block_location[2] // s * block_size[2] if nx3 > 1 else 0

                # Only process if this block contributes to the subsampled grid
                if (block_location[0] % s == 0 and
                    block_location[1] % s == 0 and
                    block_location[2] % s == 0):

                    iu_d = il_d + block_size[0] if nx1 > 1 else 1
                    ju_d = jl_d + block_size[1] if nx2 > 1 else 1
                    ku_d = kl_d + block_size[2] if nx3 > 1 else 1

                    il_s = max(il_d, i_min) - il_d
                    jl_s = max(jl_d, j_min) - jl_d
                    kl_s = max(kl_d, k_min) - kl_d
                    iu_s = min(iu_d, i_max) - il_d
                    ju_s = min(ju_d, j_max) - jl_d
                    ku_s = min(ku_d, k_max) - kl_d

                    if il_s < iu_s and jl_s < ju_s and kl_s < ku_s:
                        il_d = max(il_d, i_min) - i_min
                        jl_d = max(jl_d, j_min) - j_min
                        kl_d = max(kl_d, k_min) - k_min
                        iu_d = min(iu_d, i_max) - i_min
                        ju_d = min(ju_d, j_max) - j_min
                        ku_d = min(ku_d, k_max) - k_min

                        for q in quantities:
                            block_data = filedata['mb_data'][q][block_num]
                            # Subsample by taking every s-th point
                            subsampled = block_data[::s, ::s, ::s]
                            data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = \
                                subsampled[kl_s:ku_s, jl_s:ju_s, il_s:iu_s]

                        if return_levels:
                            data['Levels'][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level

            elif fast_restrict:
                # Fast restriction: average cells without volume weighting
                s = 2 ** (block_level - level)
                il_d = block_location[0] // s * block_size[0] if nx1 > 1 else 0
                jl_d = block_location[1] // s * block_size[1] if nx2 > 1 else 0
                kl_d = block_location[2] // s * block_size[2] if nx3 > 1 else 0
                iu_d = il_d + block_size[0] if nx1 > 1 else 1
                ju_d = jl_d + block_size[1] if nx2 > 1 else 1
                ku_d = kl_d + block_size[2] if nx3 > 1 else 1

                il_s = max(il_d, i_min) - il_d
                jl_s = max(jl_d, j_min) - jl_d
                kl_s = max(kl_d, k_min) - kl_d
                iu_s = min(iu_d, i_max) - il_d
                ju_s = min(ju_d, j_max) - jl_d
                ku_s = min(ku_d, k_max) - kl_d

                if il_s < iu_s and jl_s < ju_s and kl_s < ku_s:
                    il_d = max(il_d, i_min) - i_min
                    jl_d = max(jl_d, j_min) - j_min
                    kl_d = max(kl_d, k_min) - k_min
                    iu_d = min(iu_d, i_max) - i_min
                    ju_d = min(ju_d, j_max) - j_min
                    ku_d = min(ku_d, k_max) - k_min

                    for q in quantities:
                        block_data = filedata['mb_data'][q][block_num]
                        # Reshape and average
                        nk, nj, ni = block_data.shape
                        restricted = block_data.reshape(
                            nk // s, s, nj // s, s, ni // s, s
                        ).mean(axis=(1, 3, 5))
                        data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] += \
                            restricted[kl_s:ku_s, jl_s:ju_s, il_s:iu_s]

                    if return_levels:
                        data['Levels'][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level

            else:
                # Volume-weighted restriction (most accurate)
                s = 2 ** (block_level - level)
                loc1_coarse = block_location[0] // s
                loc2_coarse = block_location[1] // s
                loc3_coarse = block_location[2] // s

                il_d = loc1_coarse * block_size[0] if nx1 > 1 else 0
                jl_d = loc2_coarse * block_size[1] if nx2 > 1 else 0
                kl_d = loc3_coarse * block_size[2] if nx3 > 1 else 0
                iu_d = il_d + block_size[0] if nx1 > 1 else 1
                ju_d = jl_d + block_size[1] if nx2 > 1 else 1
                ku_d = kl_d + block_size[2] if nx3 > 1 else 1

                # Mark this coarse block as having restricted data
                if lx1 > 0 and lx2 > 0 and lx3 > 0:
                    restricted_data[loc3_coarse, loc2_coarse, loc1_coarse] = True

                # Compute fine cell indices within this coarse block
                sub_loc1 = block_location[0] % s
                sub_loc2 = block_location[1] % s
                sub_loc3 = block_location[2] % s

                for q in quantities:
                    block_data = filedata['mb_data'][q][block_num]
                    nk_fine, nj_fine, ni_fine = block_data.shape

                    # For each fine cell, add volume-weighted contribution
                    for kf in range(nk_fine):
                        kc = kl_d + (sub_loc3 * nk_fine + kf) // s
                        if kc < k_min or kc >= k_max:
                            continue
                        kc_out = kc - k_min

                        for jf in range(nj_fine):
                            jc = jl_d + (sub_loc2 * nj_fine + jf) // s
                            if jc < j_min or jc >= j_max:
                                continue
                            jc_out = jc - j_min

                            for if_ in range(ni_fine):
                                ic = il_d + (sub_loc1 * ni_fine + if_) // s
                                if ic < i_min or ic >= i_max:
                                    continue
                                ic_out = ic - i_min

                                # Get fine cell bounds
                                dx1 = (data['x1f'][ic + 1] - data['x1f'][ic]) / s
                                dx2 = (data['x2f'][jc + 1] - data['x2f'][jc]) / s
                                dx3 = (data['x3f'][kc + 1] - data['x3f'][kc]) / s

                                x1m = data['x1f'][ic] + (sub_loc1 * ni_fine + if_) % s * dx1
                                x1p = x1m + dx1
                                x2m = data['x2f'][jc] + (sub_loc2 * nj_fine + jf) % s * dx2
                                x2p = x2m + dx2
                                x3m = data['x3f'][kc] + (sub_loc3 * nk_fine + kf) % s * dx3
                                x3p = x3m + dx3

                                vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                data[q][kc_out, jc_out, ic_out] += \
                                    block_data[kf, jf, if_] * vol

                if return_levels:
                    for kc in range(kl_d, iu_d):
                        if kc < k_min or kc >= k_max:
                            continue
                        for jc in range(jl_d, ju_d):
                            if jc < j_min or jc >= j_max:
                                continue
                            for ic in range(il_d, iu_d):
                                if ic < i_min or ic >= i_max:
                                    continue
                                data['Levels'][kc - k_min, jc - j_min, ic - i_min] = level

    # Finalize volume-weighted restriction by dividing out volumes
    if level < max_level and not subsample and not fast_restrict:
        for loc3 in range(lx3):
            for loc2 in range(lx2):
                for loc1 in range(lx1):
                    if restricted_data[loc3, loc2, loc1]:
                        il = loc1 * block_size[0]
                        jl = loc2 * block_size[1]
                        kl = loc3 * block_size[2]
                        iu = il + block_size[0]
                        ju = jl + block_size[1]
                        ku = kl + block_size[2]
                        il = max(il, i_min) - i_min
                        jl = max(jl, j_min) - j_min
                        kl = max(kl, k_min) - k_min
                        iu = min(iu, i_max) - i_min
                        ju = min(ju, j_max) - j_min
                        ku = min(ku, k_max) - k_min

                        # Vectorized volume calculation
                        x1m = data['x1f'][il:iu][None, None, :]
                        x1p = data['x1f'][il + 1:iu + 1][None, None, :]
                        x2m = data['x2f'][jl:ju][None, :, None]
                        x2p = data['x2f'][jl + 1:ju + 1][None, :, None]
                        x3m = data['x3f'][kl:ku][:, None, None]
                        x3p = data['x3f'][kl + 1:ku + 1][:, None, None]

                        vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)

                        for q in quantities:
                            data[q][kl:ku, jl:ju, il:iu] /= vol

    # Add metadata
    data['Time'] = filedata['time']
    data['NumCycles'] = filedata['cycle']
    data['MaxLevel'] = max_level

    return data


# --------------------------------------------------------------------------------------
# Public athdf conversion functions -----------------------------------------------------
# --------------------------------------------------------------------------------------


def read_binary_as_athdf(
    filename: str,
    raw: bool = False,
    data: Optional[Dict[str, Any]] = None,
    quantities: Optional[Sequence[str]] = None,
    dtype: Optional[np.dtype] = None,
    level: Optional[int] = None,
    return_levels: bool = False,
    subsample: bool = False,
    fast_restrict: bool = False,
    x1_min: Optional[float] = None,
    x1_max: Optional[float] = None,
    x2_min: Optional[float] = None,
    x2_max: Optional[float] = None,
    x3_min: Optional[float] = None,
    x3_max: Optional[float] = None,
    vol_func: Optional[Callable] = None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1: Optional[Callable] = None,
    center_func_2: Optional[Callable] = None,
    center_func_3: Optional[Callable] = None,
    num_ghost: int = 0,
) -> Dict[str, Any]:
    """Read a bin file and organize data similar to athdf format."""
    filedata = read_binary(filename)
    if raw:
        return filedata
    return _organize_as_athdf(
        filedata, data=data, quantities=quantities, dtype=dtype, level=level,
        return_levels=return_levels, subsample=subsample, fast_restrict=fast_restrict,
        x1_min=x1_min, x1_max=x1_max, x2_min=x2_min, x2_max=x2_max,
        x3_min=x3_min, x3_max=x3_max, vol_func=vol_func,
        center_func_1=center_func_1, center_func_2=center_func_2,
        center_func_3=center_func_3, num_ghost=num_ghost,
        use_output_block_size=True,
    )


def read_all_ranks_binary_as_athdf(
    rank0_filename: str,
    raw: bool = False,
    data: Optional[Dict[str, Any]] = None,
    quantities: Optional[Sequence[str]] = None,
    dtype: Optional[np.dtype] = None,
    level: Optional[int] = None,
    return_levels: bool = False,
    subsample: bool = False,
    fast_restrict: bool = False,
    x1_min: Optional[float] = None,
    x1_max: Optional[float] = None,
    x2_min: Optional[float] = None,
    x2_max: Optional[float] = None,
    x3_min: Optional[float] = None,
    x3_max: Optional[float] = None,
    vol_func: Optional[Callable] = None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1: Optional[Callable] = None,
    center_func_2: Optional[Callable] = None,
    center_func_3: Optional[Callable] = None,
    num_ghost: int = 0,
) -> Dict[str, Any]:
    """Read binary files from all ranks and organize as athdf format."""
    filedata = read_all_ranks_binary(rank0_filename)
    if raw:
        return filedata
    return _organize_as_athdf(
        filedata, data=data, quantities=quantities, dtype=dtype, level=level,
        return_levels=return_levels, subsample=subsample, fast_restrict=fast_restrict,
        x1_min=x1_min, x1_max=x1_max, x2_min=x2_min, x2_max=x2_max,
        x3_min=x3_min, x3_max=x3_max, vol_func=vol_func,
        center_func_1=center_func_1, center_func_2=center_func_2,
        center_func_3=center_func_3, num_ghost=num_ghost,
        use_output_block_size=True,
    )


def read_coarsened_binary_as_athdf(
    filename: str,
    raw: bool = False,
    data: Optional[Dict[str, Any]] = None,
    quantities: Optional[Sequence[str]] = None,
    dtype: Optional[np.dtype] = None,
    level: Optional[int] = None,
    return_levels: bool = False,
    subsample: bool = False,
    fast_restrict: bool = False,
    x1_min: Optional[float] = None,
    x1_max: Optional[float] = None,
    x2_min: Optional[float] = None,
    x2_max: Optional[float] = None,
    x3_min: Optional[float] = None,
    x3_max: Optional[float] = None,
    vol_func: Optional[Callable] = None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1: Optional[Callable] = None,
    center_func_2: Optional[Callable] = None,
    center_func_3: Optional[Callable] = None,
    num_ghost: int = 0,
) -> Dict[str, Any]:
    """Read a coarsened bin file and organize data similar to athdf format."""
    filedata = read_coarsened_binary(filename)
    if raw:
        return filedata
    return _organize_as_athdf(
        filedata, data=data, quantities=quantities, dtype=dtype, level=level,
        return_levels=return_levels, subsample=subsample, fast_restrict=fast_restrict,
        x1_min=x1_min, x1_max=x1_max, x2_min=x2_min, x2_max=x2_max,
        x3_min=x3_min, x3_max=x3_max, vol_func=vol_func,
        center_func_1=center_func_1, center_func_2=center_func_2,
        center_func_3=center_func_3, num_ghost=num_ghost,
        use_output_block_size=False,
    )


def read_all_ranks_coarsened_binary_as_athdf(
    rank0_filename: str,
    raw: bool = False,
    data: Optional[Dict[str, Any]] = None,
    quantities: Optional[Sequence[str]] = None,
    dtype: Optional[np.dtype] = None,
    level: Optional[int] = None,
    return_levels: bool = False,
    subsample: bool = False,
    fast_restrict: bool = False,
    x1_min: Optional[float] = None,
    x1_max: Optional[float] = None,
    x2_min: Optional[float] = None,
    x2_max: Optional[float] = None,
    x3_min: Optional[float] = None,
    x3_max: Optional[float] = None,
    vol_func: Optional[Callable] = None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1: Optional[Callable] = None,
    center_func_2: Optional[Callable] = None,
    center_func_3: Optional[Callable] = None,
    num_ghost: int = 0,
) -> Dict[str, Any]:
    """Read coarsened binary files from all ranks and organize as athdf format."""
    filedata = read_all_ranks_coarsened_binary(rank0_filename)
    if raw:
        return filedata
    return _organize_as_athdf(
        filedata, data=data, quantities=quantities, dtype=dtype, level=level,
        return_levels=return_levels, subsample=subsample, fast_restrict=fast_restrict,
        x1_min=x1_min, x1_max=x1_max, x2_min=x2_min, x2_max=x2_max,
        x3_min=x3_min, x3_max=x3_max, vol_func=vol_func,
        center_func_1=center_func_1, center_func_2=center_func_2,
        center_func_3=center_func_3, num_ghost=num_ghost,
        use_output_block_size=False,
    )


def read_single_rank_binary_as_athdf(
    filename: str,
    meshblock_index_in_file: int,
    raw: bool = False,
    data: Optional[Dict[str, Any]] = None,
    quantities: Optional[Sequence[str]] = None,
    dtype: Optional[np.dtype] = None,
    return_levels: bool = False,
    x1_min: Optional[float] = None,
    x1_max: Optional[float] = None,
    x2_min: Optional[float] = None,
    x2_max: Optional[float] = None,
    x3_min: Optional[float] = None,
    x3_max: Optional[float] = None,
    vol_func: Optional[Callable] = None,
    center_func_1: Optional[Callable] = None,
    center_func_2: Optional[Callable] = None,
    center_func_3: Optional[Callable] = None,
) -> Dict[str, Any]:
    """Read a single meshblock from a rank binary file as athdf format."""
    filedata = read_binary(filename)
    if raw:
        return filedata

    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    block_size = [filedata['nx1_mb'], filedata['nx2_mb'], filedata['nx3_mb']]
    if dtype is None:
        dtype = np.float32

    if vol_func is None:
        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)
    if center_func_1 is None:
        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)
    if center_func_2 is None:
        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)
    if center_func_3 is None:
        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    center_funcs = [center_func_1, center_func_2, center_func_3]

    # Use meshblock geometry for local coordinates
    for d in range(1, 4):
        xf = f'x{d}f'
        xv = f'x{d}v'
        nx = block_size[d - 1]
        xmin = filedata['mb_geometry'][meshblock_index_in_file, (d - 1) * 2]
        xmax = filedata['mb_geometry'][meshblock_index_in_file, (d - 1) * 2 + 1]
        data[xf] = np.linspace(xmin, xmax, nx + 1, dtype=dtype)
        data[xv] = np.empty(nx, dtype=dtype)
        for i in range(nx):
            data[xv][i] = center_funcs[d - 1](data[xf][i], data[xf][i + 1])

    if quantities is None:
        quantities = filedata['var_names']

    # Account for selection
    i_min, i_max = 0, block_size[0]
    j_min, j_max = 0, block_size[1]
    k_min, k_max = 0, block_size[2]
    if x1_min is not None:
        i_min = max(i_min, np.searchsorted(data['x1f'], x1_min))
    if x1_max is not None:
        i_max = min(i_max, np.searchsorted(data['x1f'], x1_max))
    if x2_min is not None:
        j_min = max(j_min, np.searchsorted(data['x2f'], x2_min))
    if x2_max is not None:
        j_max = min(j_max, np.searchsorted(data['x2f'], x2_max))
    if x3_min is not None:
        k_min = max(k_min, np.searchsorted(data['x3f'], x3_min))
    if x3_max is not None:
        k_max = min(k_max, np.searchsorted(data['x3f'], x3_max))

    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data['Levels'] = np.empty(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)

    for q in quantities:
        block_data = filedata['mb_data'][q][meshblock_index_in_file]
        data[q] = block_data[k_min:k_max, j_min:j_max, i_min:i_max]

    if return_levels:
        data['Levels'].fill(filedata['mb_logical'][meshblock_index_in_file, 3])

    data['Time'] = filedata['time']
    data['NumCycles'] = filedata['cycle']
    data['MaxLevel'] = filedata['mb_logical'][meshblock_index_in_file, 3]

    return data


def read_rank_binary_as_athdf(
    filename: str,
    raw: bool = False,
    data: Optional[Dict[str, Any]] = None,
    quantities: Optional[Sequence[str]] = None,
    dtype: Optional[np.dtype] = None,
    return_levels: bool = False,
) -> Dict[str, Any]:
    """Read a single rank binary file and combine all meshblocks into a unified grid.

    Returns an athdf-like dictionary covering the full extent of the rank's data.
    """
    filedata = read_binary(filename)
    if raw:
        return filedata

    Nx1 = filedata['Nx1']
    Nx2 = filedata['Nx2']
    Nx3 = filedata['Nx3']
    x1min = filedata['x1min']
    x1max = filedata['x1max']
    x2min = filedata['x2min']
    x2max = filedata['x2max']
    x3min = filedata['x3min']
    x3max = filedata['x3max']
    n_mbs = filedata['n_mbs']

    if dtype is None:
        dtype = np.float32

    def center_func(xm, xp):
        return 0.5 * (xm + xp)

    x1f = np.linspace(x1min, x1max, Nx1 + 1, dtype=dtype)
    x1v = np.array([center_func(xm, xp) for xm, xp in zip(x1f[:-1], x1f[1:])], dtype=dtype)
    x2f = np.linspace(x2min, x2max, Nx2 + 1, dtype=dtype)
    x2v = np.array([center_func(xm, xp) for xm, xp in zip(x2f[:-1], x2f[1:])], dtype=dtype)
    x3f = np.linspace(x3min, x3max, Nx3 + 1, dtype=dtype)
    x3v = np.array([center_func(xm, xp) for xm, xp in zip(x3f[:-1], x3f[1:])], dtype=dtype)

    if quantities is None:
        quantities = filedata['var_names']

    data = {}
    for q in quantities:
        data[q] = np.zeros((Nx3, Nx2, Nx1), dtype=dtype)
    if return_levels:
        data['Levels'] = np.zeros((Nx3, Nx2, Nx1), dtype=np.int32)

    for mb in range(n_mbs):
        mb_index = filedata['mb_index'][mb]
        is_, ie, js, je, ks, ke = mb_index
        is_, ie = max(is_, 0), min(ie, Nx1 - 1)
        js, je = max(js, 0), min(je, Nx2 - 1)
        ks, ke = max(ks, 0), min(ke, Nx3 - 1)
        i_slice = slice(is_, ie + 1)
        j_slice = slice(js, je + 1)
        k_slice = slice(ks, ke + 1)
        for q in quantities:
            block_data = filedata['mb_data'][q][mb]
            data[q][k_slice, j_slice, i_slice] = \
                block_data[:ke - ks + 1, :je - js + 1, :ie - is_ + 1]
        if return_levels:
            level = filedata['mb_logical'][mb][3]
            data['Levels'][k_slice, j_slice, i_slice] = level

    data['x1f'] = x1f
    data['x1v'] = x1v
    data['x2f'] = x2f
    data['x2v'] = x2v
    data['x3f'] = x3f
    data['x3v'] = x3v
    data['Time'] = filedata['time']
    data['NumCycles'] = filedata['cycle']
    data['MaxLevel'] = int(np.max(filedata['mb_logical'][:, 3]))
    return data


# --------------------------------------------------------------------------------------
# XDMF and HDF5 output ------------------------------------------------------------------
# --------------------------------------------------------------------------------------


def write_xdmf_for(
    xdmfname: str,
    dumpname: str,
    fdata: Dict[str, Any],
    mode: str = "auto",
) -> None:
    """Write an XDMF file for a fluid snapshot file.

    Parameters
    ----------
    xdmfname
        Name of XDMF file to write.
    dumpname
        Location of fluid data file relative to xdmfname directory.
    fdata
        Dictionary of fluid file data, e.g., from read_binary().
    mode
        Force xdmf format (auto sets by extension).
    """
    with open(xdmfname, "w") as fp:

        def write_meshblock(fp, mb, nx1, nx2, nx3, nmb, dumpname, vars_no_b, vars_w_b):
            fp.write(f'  <Grid Name="MeshBlock{mb}" GridType="Uniform">\n')
            fp.write('   <Topology TopologyType="3DRectMesh"')
            fp.write(f' NumberOfElements="{nx3+1} {nx2+1} {nx1+1}"/>\n')
            fp.write('   <Geometry GeometryType="VXVYVZ">\n')
            fp.write(
                f'    <DataItem ItemType="HyperSlab" Dimensions="{nx1+1}">\n'
                f'     <DataItem Dimensions="3 2" NumberType="Int"> '
                f'{mb} 0 1 1 1 {nx1+1} </DataItem>\n'
                f'     <DataItem Dimensions="{nmb} {nx1+1}" Format="HDF"> '
                f'{dumpname}:/x1f </DataItem>\n'
                f'    </DataItem>\n'
                f'    <DataItem ItemType="HyperSlab" Dimensions="{nx2+1}">\n'
                f'     <DataItem Dimensions="3 2" NumberType="Int"> '
                f'{mb} 0 1 1 1 {nx2+1} </DataItem>\n'
                f'     <DataItem Dimensions="{nmb} {nx2+1}" Format="HDF"> '
                f'{dumpname}:/x2f </DataItem>\n'
                f'    </DataItem>\n'
                f'    <DataItem ItemType="HyperSlab" Dimensions="{nx3+1}">\n'
                f'     <DataItem Dimensions="3 2" NumberType="Int"> '
                f'{mb} 0 1 1 1 {nx3+1} </DataItem>\n'
                f'     <DataItem Dimensions="{nmb} {nx3+1}" Format="HDF"> '
                f'{dumpname}:/x3f </DataItem>\n'
                f'    </DataItem>\n'
                f'   </Geometry>\n'
            )

            nvar_no_b = len(vars_no_b)
            for vi, var_name in enumerate(vars_no_b):
                fp.write(
                    f'   <Attribute Name="{var_name}" Center="Cell">\n'
                    f'    <DataItem ItemType="HyperSlab" Dimensions="{nx3} {nx2} {nx1}">\n'
                    f'     <DataItem Dimensions="3 5" NumberType="Int">\n'
                    f'      {vi} {mb} 0 0 0 1 1 1 1 1 1 1 {nx3} {nx2} {nx1}\n'
                    f'     </DataItem>\n'
                    f'     <DataItem Dimensions="{nvar_no_b} {nmb} {nx3} {nx2} {nx1}" '
                    f'Format="HDF">\n'
                    f'      {dumpname}:/uov\n'
                    f'     </DataItem>\n'
                    f'    </DataItem>\n'
                    f'   </Attribute>\n'
                )

            if vars_w_b:
                nvar_w_b = len(vars_w_b)
                for vi, var_name in enumerate(vars_w_b):
                    fp.write(
                        f'   <Attribute Name="{var_name}" Center="Cell">\n'
                        f'    <DataItem ItemType="HyperSlab" Dimensions="{nx3} {nx2} {nx1}">\n'
                        f'     <DataItem Dimensions="3 5" NumberType="Int">\n'
                        f'      {vi} {mb} 0 0 0 1 1 1 1 1 1 1 {nx3} {nx2} {nx1}\n'
                        f'     </DataItem>\n'
                        f'     <DataItem Dimensions="{nvar_w_b} {nmb} {nx3} {nx2} {nx1}" '
                        f'Format="HDF">\n'
                        f'      {dumpname}:/B\n'
                        f'     </DataItem>\n'
                        f'    </DataItem>\n'
                        f'   </Attribute>\n'
                    )

            fp.write('  </Grid>\n')

        # File header
        fp.write('<?xml version="1.0" ?>\n')
        fp.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        fp.write('<Xdmf Version="2.0">\n')
        fp.write('<Information Name="TimeVaryingMetaData" Value="True"/>\n')
        fp.write('<Domain>\n')
        fp.write('<Grid Name="Mesh" GridType="Collection">\n')
        fp.write(f' <Time Value="{fdata["time"]}"/>\n')

        vars_without_b = [v for v in fdata["var_names"] if "bcc" not in v]
        vars_only_b = [v for v in fdata["var_names"] if v not in vars_without_b]

        nx1 = fdata["nx1_out_mb"]
        nx2 = fdata["nx2_out_mb"]
        nx3 = fdata["nx3_out_mb"]
        nmb = fdata["n_mbs"]

        for mb in range(nmb):
            write_meshblock(
                fp, mb, nx1, nx2, nx3, nmb, dumpname, vars_without_b, vars_only_b
            )

        fp.write('</Grid>\n')
        fp.write('</Domain>\n')
        fp.write('</Xdmf>\n')


def write_athdf(athdf_fname: str, filedata: Dict[str, Any]) -> None:
    """Write binary data to an athdf (HDF5) file.

    Parameters
    ----------
    athdf_fname
        Output filename (should end in .athdf).
    filedata
        Dictionary from read_binary() or similar.
    """
    n_mbs = filedata['n_mbs']
    nx1 = filedata['nx1_out_mb']
    nx2 = filedata['nx2_out_mb']
    nx3 = filedata['nx3_out_mb']
    var_names = filedata['var_names']
    nvars = len(var_names)

    # Separate B-field variables from others
    vars_without_b = [v for v in var_names if "bcc" not in v]
    vars_only_b = [v for v in var_names if v not in vars_without_b]

    with h5py.File(athdf_fname, 'w') as f:
        # Write metadata
        f.attrs['Time'] = filedata['time']
        f.attrs['NumCycles'] = filedata['cycle']
        f.attrs['NumMeshBlocks'] = n_mbs
        f.attrs['MaxLevel'] = int(np.max(filedata['mb_logical'][:, 3]))
        f.attrs['MeshBlockSize'] = np.array([nx1, nx2, nx3], dtype=np.int32)
        f.attrs['RootGridSize'] = np.array(
            [filedata['Nx1'], filedata['Nx2'], filedata['Nx3']], dtype=np.int32
        )
        f.attrs['RootGridX1'] = np.array(
            [filedata['x1min'], filedata['x1max']], dtype=np.float64
        )
        f.attrs['RootGridX2'] = np.array(
            [filedata['x2min'], filedata['x2max']], dtype=np.float64
        )
        f.attrs['RootGridX3'] = np.array(
            [filedata['x3min'], filedata['x3max']], dtype=np.float64
        )
        f.attrs['NumVariables'] = np.array([len(vars_without_b), len(vars_only_b)],
                                           dtype=np.int32)
        f.attrs['VariableNames'] = np.array(var_names, dtype='S')
        f.attrs['DatasetNames'] = np.array(['uov', 'B'], dtype='S')

        # Write coordinate arrays per meshblock
        x1f = np.zeros((n_mbs, nx1 + 1), dtype=np.float64)
        x2f = np.zeros((n_mbs, nx2 + 1), dtype=np.float64)
        x3f = np.zeros((n_mbs, nx3 + 1), dtype=np.float64)

        for mb in range(n_mbs):
            geom = filedata['mb_geometry'][mb]
            x1f[mb] = np.linspace(geom[0], geom[1], nx1 + 1)
            x2f[mb] = np.linspace(geom[2], geom[3], nx2 + 1)
            x3f[mb] = np.linspace(geom[4], geom[5], nx3 + 1)

        f.create_dataset('x1f', data=x1f)
        f.create_dataset('x2f', data=x2f)
        f.create_dataset('x3f', data=x3f)

        # Write logical locations and levels
        f.create_dataset('LogicalLocations', data=filedata['mb_logical'][:, :3])
        f.create_dataset('Levels', data=filedata['mb_logical'][:, 3])

        # Write variable data
        if vars_without_b:
            uov_data = np.zeros((len(vars_without_b), n_mbs, nx3, nx2, nx1),
                                dtype=np.float32)
            for vi, vname in enumerate(vars_without_b):
                for mb in range(n_mbs):
                    uov_data[vi, mb] = filedata['mb_data'][vname][mb]
            f.create_dataset('uov', data=uov_data)

        if vars_only_b:
            b_data = np.zeros((len(vars_only_b), n_mbs, nx3, nx2, nx1),
                              dtype=np.float32)
            for vi, vname in enumerate(vars_only_b):
                for mb in range(n_mbs):
                    b_data[vi, mb] = filedata['mb_data'][vname][mb]
            f.create_dataset('B', data=b_data)


def convert_file(binary_fname: str) -> None:
    """Convert a single binary file to athdf format with XDMF.

    Parameters
    ----------
    binary_fname
        Filename of bin file to convert.

    Creates "binary_data.bin" -> "binary_data.athdf" and "binary_data.athdf.xdmf"
    """
    athdf_fname = binary_fname.replace(".bin", "") + ".athdf"
    xdmf_fname = athdf_fname + ".xdmf"
    filedata = read_binary(binary_fname)
    write_athdf(athdf_fname, filedata)
    write_xdmf_for(xdmf_fname, os.path.basename(athdf_fname), filedata)


# --------------------------------------------------------------------------------------
# Athinput reader -----------------------------------------------------------------------
# --------------------------------------------------------------------------------------


def athinput(filename: str) -> Dict[str, Dict[str, Any]]:
    """Read athinput file and return a dictionary of dictionaries."""
    with open(filename, 'r') as f:
        lines = filter(None, [i.split('#')[0].strip() for i in f.readlines()])
    data: Dict[str, Dict[str, Any]] = {}
    blocks = ('\n'.join(lines)).split('<')[1:]

    def typecast(x: str) -> Any:
        if '_' in x:
            return x
        try:
            return int(x)
        except ValueError:
            pass
        try:
            return float(x)
        except ValueError:
            pass
        try:
            return complex(x)
        except ValueError:
            pass
        return x

    def parse_line(line: str) -> List[Any]:
        out = [i.strip() for i in line.split('=')]
        out[1] = '='.join(out[1:])
        out[1] = typecast(out[1])
        return out[:2]

    for block in blocks:
        info = list(filter(None, block.split('\n')))
        key = info.pop(0)[:-1]  # last character is '>'
        data[key] = dict(map(parse_line, info))
    return data


# --------------------------------------------------------------------------------------
# Public API ----------------------------------------------------------------------------
# --------------------------------------------------------------------------------------

__all__ = [
    "read_binary",
    "read_coarsened_binary",
    "read_all_ranks_binary",
    "read_all_ranks_coarsened_binary",
    "read_binary_as_athdf",
    "read_all_ranks_binary_as_athdf",
    "read_all_ranks_coarsened_binary_as_athdf",
    "read_single_rank_binary_as_athdf",
    "read_coarsened_binary_as_athdf",
    "write_xdmf_for",
    "write_athdf",
    "convert_file",
    "read_rank_binary_as_athdf",
    "athinput",
]
