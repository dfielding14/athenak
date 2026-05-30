"""
Functions to:
  (1) convert bin --> Python dictionary
  (2) convert Python dictionary --> athdf(xdmf) files

This module contains a collection of helper functions for readng and
writing athena file data formats. More information is provided in the
function docstrings.

----

In order to translate a binary file into athdf and corresponding xdmf
files, you could do the following:

  import bin_convert
  import os

  binary_fname = "path/to/file.bin"
  athdf_fname = binary_fname.replace(".bin", ".athdf")
  xdmf_fname = athdf_fname + ".xdmf"
  filedata = bin_convert.read_binary(binary_fname)
  bin_convert.write_athdf(athdf_fname, filedata)
  bin_convert.write_xdmf_for(xdmf_fname, os.path.basename(athdf_fname), filedata)

Notice that write_xdmf_for(...) function expects the relative path to
the athdf file from the xdmf, so please be aware of this requirement.

----

The read_*(...) functions return a filedata dictionary-like object with

    filedata['header'] = array of strings
        ordered array of header, including all the header information
    filedata['time'] = float
        time from input file
    filedata['cycle'] = int
        cycle from input file
    filedata['var_names'] = array of strings
        ordered array of variable names, like ['dens', 'eint', ...]
    filedata['n_mbs'] = int
        total number of meshblocks in the file
    filedata['nx1_mb'] = int
        number of cells in x1 direction in MeshBlock
    filedata['nx2_mb'] = int
        number of cells in x2 direction in MeshBlock
    filedata['nx3_mb'] = int
        number of cells in x3 direction in MeshBlock
    filedata['nx1_out_mb'] = int
        number of output cells in x1 direction in MeshBlock (useful for slicing)
    filedata['nx2_out_mb'] = int
        number of output cells in x2 direction in MeshBlock (useful for slicing)
    filedata['nx3_out_mb'] = int
        number of output cells in x3 direction in MeshBlock (useful for slicing)
    filedata['Nx1'] = int
        total number of cell in x1 direction in root grid
    filedata['Nx2'] = int
        total number of cell in x2 direction in root grid
    filedata['Nx3'] = int
        total number of cell in x3 direction in root grid
    filedata['x1min'] = float
        coordinate minimum of root grid in x1 direction
    filedata['x1max'] = float
        coordinate maximum of root grid in x1 direction
    filedata['x2min'] = float
        coordinate minimum of root grid in x2 direction
    filedata['x2max'] = float
        coordinate maximum of root grid in x2 direction
    filedata['x3min'] = float
        coordinate minimum of root grid in x3 direction
    filedata['x3max'] = float
        coordinate maximum of root grid in x3 direction
    filedata['nvars'] = int
        number of output variables (including magnetic field if it exists)
    filedata['mb_index'] = array with shape [n_mbs, 6]
        is,ie,js,je,ks,ke range for output MeshBlock indexing (useful for slicing)
    filedata['mb_logical'] = array with shape [n_mbs, 4]
        i,j,k,level coordinates for each MeshBlock
    filedata['mb_geometry'] = array with shape [n_mbs, 6]
        x1i,x2i,x3i,dx1,dx2,dx3 including cell-centered location of left-most
        cell and offsets between cells
    filedata['mb_data'] = dict of arrays with shape [n_mbs, nx3, nx2, nx1]
        {'var1':var1_array, 'var2':var2_array, ...} dictionary of fluid data arrays
        for each variable in var_names
"""

import numpy as np
import os
import h5py
import glob
from numbers import Integral
import re


_SHARD_DIRECTORY_RE = re.compile(r"^(rank|node)_([0-9]+)$")


def _partition_info(filename):
    """Return a strict ``(kind, integer ID)`` pair for a shard path."""
    shard_name = os.path.basename(os.path.dirname(os.path.abspath(filename)))
    match = _SHARD_DIRECTORY_RE.fullmatch(shard_name)
    if match is not None:
        return match.group(1), int(match.group(2))
    if shard_name.startswith(("rank_", "node_")):
        raise ValueError(f"invalid binary shard directory {shard_name!r}")
    return "shared", None


def _is_partitioned_path(filename):
    """Return whether *filename* is under a rank_*/node_* shard directory."""
    return _partition_info(filename)[0] != "shared"


def _glob_partition_files(shard_filename):
    """Return all sibling files belonging to a rank- or node-sharded output."""
    shard_filename = os.path.abspath(shard_filename)
    shard_dir = os.path.dirname(shard_filename)
    parent_dir = os.path.dirname(shard_dir)
    base_name = os.path.basename(shard_filename)
    shard_kind, _ = _partition_info(shard_filename)
    if shard_kind == "rank":
        pattern = os.path.join(parent_dir, "rank_*", base_name)
    elif shard_kind == "node":
        pattern = os.path.join(parent_dir, "node_*", base_name)
    else:
        return [shard_filename]

    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"no binary shard files found for {shard_filename!r} "
            f"(pattern {pattern!r})"
        )
    for candidate in files:
        candidate_kind, _ = _partition_info(candidate)
        if candidate_kind != shard_kind:
            raise ValueError(
                f"binary shard {candidate!r} does not match {shard_kind!r} inventory"
            )
    return files


def _optional_pheader_int(pheader, key, filename, family):
    """Parse one additive binary preheader integer without requiring it in old files."""
    if key not in pheader:
        return None
    try:
        return int(pheader[key])
    except ValueError as exc:
        raise ValueError(
            f"{family} file {filename!r} has invalid {key!r} metadata"
        ) from exc


def _binary_partition_metadata(filename, pheader, family):
    """Validate optional additive node-shard metadata while preserving old files."""
    kind, shard_id = _partition_info(filename)
    distribution = pheader.get("distribution")
    if distribution is not None:
        if distribution not in ("shared", "rank", "node"):
            raise ValueError(
                f"{family} file {filename!r} has invalid distribution {distribution!r}"
            )
        if distribution != kind:
            raise ValueError(
                f"{family} file {filename!r} declares distribution={distribution!r} "
                f"but resides in a {kind!r} layout"
            )
    node = _optional_pheader_int(pheader, "node", filename, family)
    number_of_nodes = _optional_pheader_int(
        pheader, "number of nodes", filename, family
    )
    number_of_meshblocks = _optional_pheader_int(
        pheader, "number of meshblocks", filename, family
    )
    if node is not None:
        if kind != "node" or node != shard_id:
            raise ValueError(
                f"{family} file {filename!r} declares node={node}, "
                f"but resides in {os.path.basename(os.path.dirname(filename))!r}"
            )
    if number_of_nodes is not None:
        if kind != "node" or number_of_nodes <= 0:
            raise ValueError(
                f"{family} file {filename!r} has invalid 'number of nodes' metadata"
            )
    if number_of_meshblocks is not None and number_of_meshblocks < 0:
        raise ValueError(
            f"{family} file {filename!r} has invalid 'number of meshblocks' metadata"
        )
    return {
        "distribution": distribution if distribution is not None else kind,
        "node": node,
        "number_of_nodes": number_of_nodes,
        "number_of_meshblocks": number_of_meshblocks,
    }


def _validate_binary_sibling_inventory(shard_files, shard_data, family):
    """Require a complete unique node-ID inventory when new count metadata exists."""
    node_counts = [item.get("number_of_nodes") for item in shard_data]
    if all(count is None for count in node_counts):
        return
    if any(count is None for count in node_counts) or len(set(node_counts)) != 1:
        raise ValueError(
            f"{family} node shard inventory has inconsistent 'number of nodes' metadata"
        )
    expected_count = node_counts[0]
    sibling_ids = []
    for path, item in zip(shard_files, shard_data):
        kind, shard_id = _partition_info(path)
        if kind != "node" or item.get("distribution") != "node":
            raise ValueError(
                f"{family} node shard inventory includes non-node file {path!r}"
            )
        if item.get("node") is None:
            raise ValueError(f"{family} node shard {path!r} is missing node metadata")
        if item["node"] != shard_id:
            raise ValueError(
                f"{family} node shard {path!r} declares node={item['node']}, "
                f"but its directory identifies {shard_id}"
            )
        sibling_ids.append(shard_id)
    if len(set(sibling_ids)) != len(sibling_ids):
        raise ValueError(f"{family} node shard inventory contains duplicate node IDs")
    if expected_count != len(sibling_ids):
        raise ValueError(
            f"{family} node shard inventory is incomplete: "
            f"expected {expected_count} shards, found {len(sibling_ids)}"
        )
    expected_ids = set(range(len(sibling_ids)))
    actual_ids = set(sibling_ids)
    if actual_ids != expected_ids:
        raise ValueError(
            f"{family} node shard inventory is incomplete: "
            f"expected IDs {sorted(expected_ids)!r}, found {sorted(actual_ids)!r}"
        )


def _validate_binary_shard(reference, candidate, path, family):
    """Validate metadata that must agree before binary shards can be combined."""
    keys = (
        "time",
        "cycle",
        "var_names",
        "nvars",
        "Nx1",
        "Nx2",
        "Nx3",
        "x1min",
        "x1max",
        "x2min",
        "x2max",
        "x3min",
        "x3max",
        "nx1_mb",
        "nx2_mb",
        "nx3_mb",
    )
    if family == "coarsened binary":
        keys += ("number_of_moments",)
    for key in keys:
        if candidate[key] != reference[key]:
            raise ValueError(
                f"{family} shard metadata mismatch for {key!r} in {path!r}: "
                f"{candidate[key]!r} != {reference[key]!r}"
            )
    if candidate["header"] != reference["header"]:
        raise ValueError(f"{family} shard metadata mismatch for 'header' in {path!r}")


def _combine_partitioned_binary(shard_filename, reader, family):
    """Combine one binary shard family, retaining valid empty sliced shards."""
    shard_files = _glob_partition_files(shard_filename)
    shard_data = [reader(path) for path in shard_files]
    reference = shard_data[0]
    for path, candidate in zip(shard_files[1:], shard_data[1:]):
        _validate_binary_shard(reference, candidate, path, family)
    _validate_binary_sibling_inventory(shard_files, shard_data, family)

    nonempty = next((item for item in shard_data if item["n_mbs"] > 0), reference)
    nonempty_shape = tuple(nonempty[f"nx{axis}_out_mb"] for axis in (1, 2, 3))
    for path, candidate in zip(shard_files, shard_data):
        if candidate["n_mbs"] == 0:
            continue
        candidate_shape = tuple(candidate[f"nx{axis}_out_mb"] for axis in (1, 2, 3))
        if candidate_shape != nonempty_shape:
            raise ValueError(
                f"{family} shard output-shape mismatch in {path!r}: "
                f"{candidate_shape!r} != {nonempty_shape!r}"
            )
    combined = nonempty.copy()
    combined["mb_index"] = []
    combined["mb_logical"] = []
    combined["mb_geometry"] = []
    combined["mb_data"] = {var: [] for var in reference["var_names"]}
    logical_owners = {}
    for path, item in zip(shard_files, shard_data):
        for logical in item["mb_logical"]:
            logical_key = tuple(int(value) for value in logical)
            if logical_key in logical_owners:
                raise ValueError(
                    f"{family} shards contain duplicate logical MeshBlock "
                    f"{logical_key!r} in {path!r} and {logical_owners[logical_key]!r}"
                )
            logical_owners[logical_key] = path
        combined["mb_index"].extend(item["mb_index"])
        combined["mb_logical"].extend(item["mb_logical"])
        combined["mb_geometry"].extend(item["mb_geometry"])
        for var in reference["var_names"]:
            combined["mb_data"][var].extend(item["mb_data"][var])

    combined["mb_index"] = np.asarray(combined["mb_index"], dtype=np.int64)
    combined["mb_logical"] = np.asarray(combined["mb_logical"], dtype=np.int32)
    combined["mb_geometry"] = np.asarray(combined["mb_geometry"])
    for var in reference["var_names"]:
        combined["mb_data"][var] = np.asarray(combined["mb_data"][var])
    combined["n_mbs"] = len(combined["mb_index"])
    meshblock_counts = [item["number_of_meshblocks"] for item in shard_data]
    combined["number_of_meshblocks"] = (
        sum(meshblock_counts)
        if all(count is not None for count in meshblock_counts)
        else None
    )
    combined["shard_files"] = shard_files
    combined["distribution"] = os.path.basename(
        os.path.dirname(os.path.abspath(shard_filename))
        ).split("_", 1)[0] if _is_partitioned_path(shard_filename) else "shared"
    return combined


def _require_meshblocks_for_athdf(filedata, filename):
    """Reject athdf-like conversion when a valid binary input has no MeshBlocks."""
    if filedata["n_mbs"] == 0:
        raise ValueError(
            f"cannot convert {filename!r} to athdf-like data: "
            "binary output contains no meshblocks"
        )


def _read_meshblocks(fp, filesize, nghost, locsizebytes, varsizebytes, var_list, family):
    """Read complete meshblock records and reject partial binary payloads."""
    dtype_loc = np.float64 if locsizebytes == 8 else np.float32
    dtype_var = np.float64 if varsizebytes == 8 else np.float32
    nvars = len(var_list)
    mb_index = []
    mb_logical = []
    mb_geometry = []
    mb_data = {var: [] for var in var_list}
    fixed_bytes = 24 + 16 + 6 * locsizebytes

    while fp.tell() < filesize:
        if filesize - fp.tell() < fixed_bytes:
            raise ValueError(
                f"truncated {family} meshblock metadata in {fp.name!r}"
            )
        index = np.frombuffer(fp.read(24), dtype=np.int32).astype(np.int64) - nghost
        logical = np.frombuffer(fp.read(16), dtype=np.int32)
        geometry = np.frombuffer(fp.read(6 * locsizebytes), dtype=dtype_loc)
        shape = (
            int(index[5] - index[4] + 1),
            int(index[3] - index[2] + 1),
            int(index[1] - index[0] + 1),
        )
        if any(length <= 0 for length in shape):
            raise ValueError(
                f"invalid {family} meshblock extent {shape!r} in {fp.name!r}"
            )
        value_count = nvars * int(np.prod(shape))
        value_bytes = value_count * varsizebytes
        raw_values = fp.read(value_bytes)
        if len(raw_values) != value_bytes:
            raise ValueError(
                f"truncated {family} meshblock values in {fp.name!r}: "
                f"expected {value_bytes} bytes, found {len(raw_values)}"
            )
        values = np.frombuffer(raw_values, dtype=dtype_var).reshape((nvars,) + shape)
        mb_index.append(index)
        mb_logical.append(logical)
        mb_geometry.append(geometry)
        for vari, var in enumerate(var_list):
            mb_data[var].append(values[vari])

    return mb_index, mb_logical, mb_geometry, mb_data


def read_binary(filename, assemble_shards=False):
    """
    Reads a bin file from filename to dictionary.

    Originally written by Lev Arzamasskiy (leva@ias.edu) on 11/15/2021
    Updated to support mesh refinement by George Wong (gnwong@ias.edu) on 01/27/2022
    Made faster by Drummond Fielding on 09/09/2024

    args:
      filename - string
          filename of bin file to read
      assemble_shards - bool, optional
          when True and filename is in a rank_* or node_* directory, discover
          sibling shards and return their combined meshblocks

    returns:
      filedata - dict
          dictionary of fluid file data
    """

    if assemble_shards and _is_partitioned_path(filename):
        return _combine_partitioned_binary(
            filename, lambda path: read_binary(path), "binary"
        )

    filedata = {}

    with open(filename, "rb") as fp:
        # load file and get size
        fp.seek(0, 2)
        filesize = fp.tell()
        fp.seek(0, 0)

        # load header information and validate file format
        code_header = fp.readline().split()
        if len(code_header) < 1:
            raise TypeError("unknown file format")
        if code_header[0] != b"Athena":
            raise TypeError(
                f"bad file format \"{code_header[0].decode('utf-8')}\" "
                + '(should be "Athena")'
            )
        version = code_header[-1].split(b"=")[-1]
        if version != b"1.1":
            raise TypeError(
                f"unsupported file format version {version.decode('utf-8')}"
            )

        pheader_count = int(fp.readline().split(b"=")[-1])
        pheader = {}
        for _ in range(pheader_count - 1):
            key, val = [x.strip() for x in fp.readline().decode("utf-8").split("=")]
            pheader[key] = val
        time = float(pheader["time"])
        cycle = int(pheader["cycle"])
        locsizebytes = int(pheader["size of location"])
        varsizebytes = int(pheader["size of variable"])
        partition_metadata = _binary_partition_metadata(filename, pheader, "binary")

        nvars = int(fp.readline().split(b"=")[-1])
        var_list = [v.decode("utf-8") for v in fp.readline().split()[1:]]
        header_size = int(fp.readline().split(b"=")[-1])
        header = [
            line.decode("utf-8").split("#")[0].strip()
            for line in fp.read(header_size).split(b"\n")
        ]
        header = [line for line in header if len(line) > 0]

        if locsizebytes not in [4, 8]:
            raise ValueError(f"unsupported location size (in bytes) {locsizebytes}")
        if varsizebytes not in [4, 8]:
            raise ValueError(f"unsupported variable size (in bytes) {varsizebytes}")

        # load grid information from header and validate
        def get_from_header(header, blockname, keyname):
            blockname = blockname.strip()
            keyname = keyname.strip()
            if not blockname.startswith("<"):
                blockname = "<" + blockname
            if blockname[-1] != ">":
                blockname += ">"
            block = "<none>"
            for line in [entry for entry in header]:
                if line.startswith("<"):
                    block = line
                    continue
                key, value = line.split("=")
                if block == blockname and key.strip() == keyname:
                    return value
            raise KeyError(f"no parameter called {blockname}/{keyname}")

        Nx1 = int(get_from_header(header, "<mesh>", "nx1"))
        Nx2 = int(get_from_header(header, "<mesh>", "nx2"))
        Nx3 = int(get_from_header(header, "<mesh>", "nx3"))
        nx1 = int(get_from_header(header, "<meshblock>", "nx1"))
        nx2 = int(get_from_header(header, "<meshblock>", "nx2"))
        nx3 = int(get_from_header(header, "<meshblock>", "nx3"))

        nghost = int(get_from_header(header, "<mesh>", "nghost"))

        x1min = float(get_from_header(header, "<mesh>", "x1min"))
        x1max = float(get_from_header(header, "<mesh>", "x1max"))
        x2min = float(get_from_header(header, "<mesh>", "x2min"))
        x2max = float(get_from_header(header, "<mesh>", "x2max"))
        x3min = float(get_from_header(header, "<mesh>", "x3min"))
        x3max = float(get_from_header(header, "<mesh>", "x3max"))

        if len(var_list) != nvars:
            raise ValueError(
                f"binary variable count mismatch in {filename!r}: "
                f"declared {nvars}, listed {len(var_list)}"
            )
        mb_index, mb_logical, mb_geometry, mb_data = _read_meshblocks(
            fp, filesize, nghost, locsizebytes, varsizebytes, var_list, "binary"
        )
        mb_count = len(mb_index)
        if (
            partition_metadata["number_of_meshblocks"] is not None
            and partition_metadata["number_of_meshblocks"] != mb_count
        ):
            raise ValueError(
                f"binary file {filename!r} declares "
                f"{partition_metadata['number_of_meshblocks']} meshblocks, "
                f"found {mb_count}"
            )

    filedata["header"] = header
    filedata["time"] = time
    filedata["cycle"] = cycle
    filedata["var_names"] = var_list

    filedata["Nx1"] = Nx1
    filedata["Nx2"] = Nx2
    filedata["Nx3"] = Nx3
    filedata["nvars"] = nvars

    filedata["x1min"] = x1min
    filedata["x1max"] = x1max
    filedata["x2min"] = x2min
    filedata["x2max"] = x2max
    filedata["x3min"] = x3min
    filedata["x3max"] = x3max

    filedata["n_mbs"] = mb_count
    filedata["nx1_mb"] = nx1
    filedata["nx2_mb"] = nx2
    filedata["nx3_mb"] = nx3
    filedata["nx1_out_mb"] = (mb_index[0][1] - mb_index[0][0]) + 1 if mb_index else 0
    filedata["nx2_out_mb"] = (mb_index[0][3] - mb_index[0][2]) + 1 if mb_index else 0
    filedata["nx3_out_mb"] = (mb_index[0][5] - mb_index[0][4]) + 1 if mb_index else 0

    filedata["mb_index"] = np.array(mb_index)
    filedata["mb_logical"] = np.array(mb_logical)
    filedata["mb_geometry"] = np.array(mb_geometry)
    filedata["mb_data"] = mb_data
    filedata.update(partition_metadata)

    return filedata


def read_coarsened_binary(filename, assemble_shards=False):
    """
    Reads a coarsened bin file from filename to dictionary.
    Originally written by Lev Arzamasskiy (leva@ias.edu) on 11/15/2021
    Updated to support mesh refinement by George Wong (gnwong@ias.edu) on 01/27/2022
    Updated to support coarsened outputs and for speed by Drummond Fielding on 09/09/2024

    args:
      filename - string
          filename of bin file to read
      assemble_shards - bool, optional
          when True and filename is in a rank_* or node_* directory, discover
          sibling shards and return their combined meshblocks

    returns:
      filedata - dict
          dictionary of fluid file data
    """

    if assemble_shards and _is_partitioned_path(filename):
        return _combine_partitioned_binary(
            filename, lambda path: read_coarsened_binary(path), "coarsened binary"
        )

    filedata = {}

    with open(filename, "rb") as fp:
        # load file and get size
        fp.seek(0, 2)
        filesize = fp.tell()
        fp.seek(0, 0)

        # load header information and validate file format
        code_header = fp.readline().split()
        if len(code_header) < 1:
            raise TypeError("unknown file format")
        if code_header[0] != b"Athena":
            raise TypeError(
                f"bad file format \"{code_header[0].decode('utf-8')}\" "
                + '(should be "Athena")'
            )
        version = code_header[-1].split(b"=")[-1]
        if version != b"1.1":
            raise TypeError(
                f"unsupported file format version {version.decode('utf-8')}"
            )

        pheader_count = int(fp.readline().split(b"=")[-1])
        pheader = {}
        for _ in range(pheader_count - 1):
            key, val = [x.strip() for x in fp.readline().decode("utf-8").split("=")]
            pheader[key] = val
        time = float(pheader["time"])
        cycle = int(pheader["cycle"])
        locsizebytes = int(pheader["size of location"])
        varsizebytes = int(pheader["size of variable"])
        coarsen_factor = int(pheader["coarsening factor"])
        partition_metadata = _binary_partition_metadata(
            filename, pheader, "coarsened binary"
        )

        nvars = int(fp.readline().split(b"=")[-1])
        var_list = [v.decode("utf-8") for v in fp.readline().split()[1:]]
        header_size = int(fp.readline().split(b"=")[-1])
        header = [
            line.decode("utf-8").split("#")[0].strip()
            for line in fp.read(header_size).split(b"\n")
        ]
        header = [line for line in header if len(line) > 0]

        if locsizebytes not in [4, 8]:
            raise ValueError(f"unsupported location size (in bytes) {locsizebytes}")
        if varsizebytes not in [4, 8]:
            raise ValueError(f"unsupported variable size (in bytes) {varsizebytes}")

        # load grid information from header and validate
        def get_from_header(header, blockname, keyname):
            blockname = blockname.strip()
            keyname = keyname.strip()
            if not blockname.startswith("<"):
                blockname = "<" + blockname
            if blockname[-1] != ">":
                blockname += ">"
            block = "<none>"
            for line in [entry for entry in header]:
                if line.startswith("<"):
                    block = line
                    continue
                key, value = line.split("=")
                if block == blockname and key.strip() == keyname:
                    return value
            raise KeyError(f"no parameter called {blockname}/{keyname}")

        Nx1 = int(get_from_header(header, "<mesh>", "nx1"))
        Nx2 = int(get_from_header(header, "<mesh>", "nx2"))
        Nx3 = int(get_from_header(header, "<mesh>", "nx3"))
        nx1 = int(get_from_header(header, "<meshblock>", "nx1"))
        nx2 = int(get_from_header(header, "<meshblock>", "nx2"))
        nx3 = int(get_from_header(header, "<meshblock>", "nx3"))

        nghost = int(get_from_header(header, "<mesh>", "nghost"))

        x1min = float(get_from_header(header, "<mesh>", "x1min"))
        x1max = float(get_from_header(header, "<mesh>", "x1max"))
        x2min = float(get_from_header(header, "<mesh>", "x2min"))
        x2max = float(get_from_header(header, "<mesh>", "x2max"))
        x3min = float(get_from_header(header, "<mesh>", "x3min"))
        x3max = float(get_from_header(header, "<mesh>", "x3max"))

        if len(var_list) != nvars:
            raise ValueError(
                f"coarsened binary variable count mismatch in {filename!r}: "
                f"declared {nvars}, listed {len(var_list)}"
            )
        mb_index, mb_logical, mb_geometry, mb_data = _read_meshblocks(
            fp,
            filesize,
            nghost,
            locsizebytes,
            varsizebytes,
            var_list,
            "coarsened binary",
        )
        mb_count = len(mb_index)
        if (
            partition_metadata["number_of_meshblocks"] is not None
            and partition_metadata["number_of_meshblocks"] != mb_count
        ):
            raise ValueError(
                f"coarsened binary file {filename!r} declares "
                f"{partition_metadata['number_of_meshblocks']} meshblocks, "
                f"found {mb_count}"
            )

    filedata["header"] = header
    filedata["time"] = time
    filedata["cycle"] = cycle
    filedata["var_names"] = var_list

    filedata["Nx1"] = Nx1 // coarsen_factor
    filedata["Nx2"] = Nx2 // coarsen_factor
    filedata["Nx3"] = Nx3 // coarsen_factor
    filedata["nvars"] = nvars
    filedata["number_of_moments"] = int(pheader["number of moments"])

    filedata["x1min"] = x1min
    filedata["x1max"] = x1max
    filedata["x2min"] = x2min
    filedata["x2max"] = x2max
    filedata["x3min"] = x3min
    filedata["x3max"] = x3max

    filedata["n_mbs"] = mb_count
    filedata["nx1_mb"] = nx1 // coarsen_factor
    filedata["nx2_mb"] = nx2 // coarsen_factor
    filedata["nx3_mb"] = nx3 // coarsen_factor
    filedata["nx1_out_mb"] = (mb_index[0][1] - mb_index[0][0]) + 1 if mb_index else 0
    filedata["nx2_out_mb"] = (mb_index[0][3] - mb_index[0][2]) + 1 if mb_index else 0
    filedata["nx3_out_mb"] = (mb_index[0][5] - mb_index[0][4]) + 1 if mb_index else 0

    filedata["mb_index"] = np.array(mb_index)
    filedata["mb_logical"] = np.array(mb_logical)
    filedata["mb_geometry"] = np.array(mb_geometry)
    filedata["mb_data"] = mb_data
    filedata.update(partition_metadata)

    return filedata


def read_all_ranks_binary(rank0_filename):
    """
    Reads binary files from all rank or node shards into a single dictionary.

    args:
      rank0_filename - string
          filename of any rank/node shard, or a shared binary file

    returns:
      combined_filedata - dict
          dictionary of combined fluid file data from all ranks
    """
    return _combine_partitioned_binary(rank0_filename, read_binary, "binary")


def read_all_ranks_coarsened_binary(rank0_filename):
    """
    Reads coarsened binary files from all rank or node shards.

    args:
      rank0_filename - string
          filename of any rank/node shard, or a shared coarsened binary file

    returns:
      combined_filedata - dict
          dictionary of combined fluid file data from all ranks
    """
    return _combine_partitioned_binary(
        rank0_filename, read_coarsened_binary, "coarsened binary"
    )


def read_binary_as_athdf(
    filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
):
    """
    Reads a bin file and organizes data similar to athdf format without writing to file.
    """
    # Step 1: Read binary data
    filedata = read_binary(filename)

    # Step 2: Organize data similar to athdf
    if raw:
        return filedata

    _require_meshblocks_for_athdf(filedata, filename)

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    max_level = max(filedata["mb_logical"][:, 3])
    if level is None:
        level = max_level
    block_size = [
        filedata["nx1_out_mb"],
        filedata["nx2_out_mb"],
        filedata["nx3_out_mb"],
    ]
    root_grid_size = [filedata["Nx1"], filedata["Nx2"], filedata["Nx3"]]
    levels = filedata["mb_logical"][:, 3]
    logical_locations = filedata["mb_logical"][:, :3]
    if dtype is None:
        dtype = np.float32

    # Calculate nx_vals
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:  # sum or slice
            other_locations = [
                location
                for location in zip(
                    levels,
                    logical_locations[:, (d + 1) % 3],
                    logical_locations[:, (d + 2) % 3],
                )
            ]
            if len(set(other_locations)) == len(other_locations):  # effective slice
                nx_vals.append(1)
            else:  # nontrivial sum
                num_blocks_this_dim = 0
                for level_this_dim, loc_this_dim in zip(
                    levels, logical_locations[:, d]
                ):
                    if level_this_dim <= level:
                        possible_max = (loc_this_dim + 1) * 2 ** (
                            level - level_this_dim
                        )
                        num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                    else:
                        possible_max = (loc_this_dim + 1) // 2 ** (
                            level_this_dim - level
                        )
                        num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                nx_vals.append(num_blocks_this_dim)
        elif block_size[d] == 1:  # singleton dimension
            nx_vals.append(1)
        else:  # normal case
            nx_vals.append(root_grid_size[d] * 2**level + 2 * num_ghost)
    nx1, nx2, nx3 = nx_vals
    lx1, lx2, lx3 = [nx // bs for nx, bs in zip(nx_vals, block_size)]

    # Set coordinate system and related functions
    # coord = "cartesian"  # Adjust based on your data
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    # Define center functions if not provided
    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    # Populate coordinate arrays
    center_funcs = [center_func_1, center_func_2, center_func_3]
    for d in range(1, 4):
        xf = f"x{d}f"
        xv = f"x{d}v"
        nx = nx_vals[d - 1]
        if nx == 1:
            xmin = filedata[f"x{d}min"]
            xmax = filedata[f"x{d}max"]
            data[xf] = np.array([xmin, xmax], dtype=dtype)
        else:
            xmin = filedata[f"x{d}min"]
            xmax = filedata[f"x{d}max"]
            data[xf] = np.linspace(xmin, xmax, nx + 1, dtype=dtype)
        data[xv] = np.empty(nx, dtype=dtype)
        for i in range(nx):
            data[xv][i] = center_funcs[d - 1](data[xf][i], data[xf][i + 1])

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max = 0, nx1
    j_min, j_max = 0, nx2
    k_min, k_max = 0, nx3
    if x1_min is not None:
        i_min = max(i_min, np.searchsorted(data["x1f"], x1_min))
    if x1_max is not None:
        i_max = min(i_max, np.searchsorted(data["x1f"], x1_max))
    if x2_min is not None:
        j_min = max(j_min, np.searchsorted(data["x2f"], x2_min))
    if x2_max is not None:
        j_max = min(j_max, np.searchsorted(data["x2f"], x2_max))
    if x3_min is not None:
        k_min = max(k_min, np.searchsorted(data["x3f"], x3_min))
    if x3_max is not None:
        k_max = min(k_max, np.searchsorted(data["x3f"], x3_max))

    # Prepare arrays for data and bookkeeping
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.empty(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)
    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Step 3: Process each block
    for block_num in range(filedata["n_mbs"]):
        block_level = levels[block_num]
        block_location = logical_locations[block_num]

        # Implement logic for prolongation, restriction, subsampling
        if block_level <= level:
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
                block_data = filedata["mb_data"][q][block_num]
                if s > 1:
                    block_data = np.repeat(
                        np.repeat(np.repeat(block_data, s, axis=2), s, axis=1),
                        s,
                        axis=0,
                    )
                data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_data[
                    kl_s:ku_s, jl_s:ju_s, il_s:iu_s
                ]
        else:
            # Implement restriction logic here (similar to athdf function)
            pass

        if return_levels:
            data["Levels"][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level

    # Step 4: Finalize data
    if level < max_level and not subsample and not fast_restrict:
        # Remove volume factors from restricted data
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
                        for k in range(kl, ku):
                            for j in range(jl, ju):
                                for i in range(il, iu):
                                    x1m, x1p = data["x1f"][i], data["x1f"][i + 1]
                                    x2m, x2p = data["x2f"][j], data["x2f"][j + 1]
                                    x3m, x3p = data["x3f"][k], data["x3f"][k + 1]
                                    vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                    for q in quantities:
                                        data[q][k, j, i] /= vol

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = max_level

    return data


def read_rank_binary_as_athdf(
    filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
):
    """Read one binary shard through the canonical logical-location mapper."""
    return read_binary_as_athdf(**locals())


def read_all_ranks_binary_as_athdf(
    rank0_filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
):
    """
    Reads a bin file and organizes data similar to athdf format without writing to file.
    """
    # Step 1: Read binary data
    filedata = read_all_ranks_binary(rank0_filename)

    # Step 2: Organize data similar to athdf
    if raw:
        return filedata

    _require_meshblocks_for_athdf(filedata, rank0_filename)

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    max_level = max(filedata["mb_logical"][:, 3])
    if level is None:
        level = max_level
    block_size = [
        filedata["nx1_out_mb"],
        filedata["nx2_out_mb"],
        filedata["nx3_out_mb"],
    ]
    root_grid_size = [filedata["Nx1"], filedata["Nx2"], filedata["Nx3"]]
    levels = filedata["mb_logical"][:, 3]
    logical_locations = filedata["mb_logical"][:, :3]
    if dtype is None:
        dtype = np.float32

    # Calculate nx_vals
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:  # sum or slice
            other_locations = [
                location
                for location in zip(
                    levels,
                    logical_locations[:, (d + 1) % 3],
                    logical_locations[:, (d + 2) % 3],
                )
            ]
            if len(set(other_locations)) == len(other_locations):  # effective slice
                nx_vals.append(1)
            else:  # nontrivial sum
                num_blocks_this_dim = 0
                for level_this_dim, loc_this_dim in zip(
                    levels, logical_locations[:, d]
                ):
                    if level_this_dim <= level:
                        possible_max = (loc_this_dim + 1) * 2 ** (
                            level - level_this_dim
                        )
                        num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                    else:
                        possible_max = (loc_this_dim + 1) // 2 ** (
                            level_this_dim - level
                        )
                        num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                nx_vals.append(num_blocks_this_dim)
        elif block_size[d] == 1:  # singleton dimension
            nx_vals.append(1)
        else:  # normal case
            nx_vals.append(root_grid_size[d] * 2**level + 2 * num_ghost)
    nx1, nx2, nx3 = nx_vals
    lx1, lx2, lx3 = [nx // bs for nx, bs in zip(nx_vals, block_size)]
    # Set coordinate system and related functions
    # coord = "cartesian"  # Adjust based on your data
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    # Define center functions if not provided
    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    # Populate coordinate arrays
    center_funcs = [center_func_1, center_func_2, center_func_3]
    for d in range(1, 4):
        xf = f"x{d}f"
        xv = f"x{d}v"
        nx = nx_vals[d - 1]
        if nx == 1:
            xmin = filedata[f"x{d}min"]
            xmax = filedata[f"x{d}max"]
            data[xf] = np.array([xmin, xmax], dtype=dtype)
        else:
            xmin = filedata[f"x{d}min"]
            xmax = filedata[f"x{d}max"]
            data[xf] = np.linspace(xmin, xmax, nx + 1, dtype=dtype)
        data[xv] = np.empty(nx, dtype=dtype)
        for i in range(nx):
            data[xv][i] = center_funcs[d - 1](data[xf][i], data[xf][i + 1])

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max = 0, nx1
    j_min, j_max = 0, nx2
    k_min, k_max = 0, nx3
    if x1_min is not None:
        i_min = max(i_min, np.searchsorted(data["x1f"], x1_min))
    if x1_max is not None:
        i_max = min(i_max, np.searchsorted(data["x1f"], x1_max))
    if x2_min is not None:
        j_min = max(j_min, np.searchsorted(data["x2f"], x2_min))
    if x2_max is not None:
        j_max = min(j_max, np.searchsorted(data["x2f"], x2_max))
    if x3_min is not None:
        k_min = max(k_min, np.searchsorted(data["x3f"], x3_min))
    if x3_max is not None:
        k_max = min(k_max, np.searchsorted(data["x3f"], x3_max))

    # Prepare arrays for data and bookkeeping
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.empty(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)
    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Step 3: Process each block
    for block_num in range(filedata["n_mbs"]):
        block_level = levels[block_num]
        block_location = logical_locations[block_num]

        # Implement logic for prolongation, restriction, subsampling
        if block_level <= level:
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
                block_data = filedata["mb_data"][q][block_num]
                if s > 1:
                    block_data = np.repeat(
                        np.repeat(np.repeat(block_data, s, axis=2), s, axis=1),
                        s,
                        axis=0,
                    )
                data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_data[
                    kl_s:ku_s, jl_s:ju_s, il_s:iu_s
                ]
        else:
            # Implement restriction logic here (similar to athdf function)
            pass

        if return_levels:
            data["Levels"][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level

    # Step 4: Finalize data
    if level < max_level and not subsample and not fast_restrict:
        # Remove volume factors from restricted data
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
                        for k in range(kl, ku):
                            for j in range(jl, ju):
                                for i in range(il, iu):
                                    x1m, x1p = data["x1f"][i], data["x1f"][i + 1]
                                    x2m, x2p = data["x2f"][j], data["x2f"][j + 1]
                                    x3m, x3p = data["x3f"][k], data["x3f"][k + 1]
                                    vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                    for q in quantities:
                                        data[q][k, j, i] /= vol

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = max_level

    return data


def read_all_ranks_coarsened_binary_as_athdf(
    rank0_filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
):
    """
    Reads a bin file and organizes data similar to athdf format without writing to file.
    """
    # Step 1: Read binary data
    filedata = read_all_ranks_coarsened_binary(rank0_filename)

    # Step 2: Organize data similar to athdf
    if raw:
        return filedata

    _require_meshblocks_for_athdf(filedata, rank0_filename)

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    max_level = max(filedata["mb_logical"][:, 3])
    if level is None:
        level = max_level
    block_size = [filedata["nx1_mb"], filedata["nx2_mb"], filedata["nx3_mb"]]
    root_grid_size = [filedata["Nx1"], filedata["Nx2"], filedata["Nx3"]]
    levels = filedata["mb_logical"][:, 3]
    logical_locations = filedata["mb_logical"][:, :3]
    if dtype is None:
        dtype = np.float32

    # Calculate nx_vals
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:
            # Implement logic for sum or slice as in athdf
            nx_vals.append(root_grid_size[d] * 2**level)
        elif block_size[d] == 1:
            nx_vals.append(1)
        else:
            nx_vals.append(root_grid_size[d] * 2**level + 2 * num_ghost)
    nx1, nx2, nx3 = nx_vals
    lx1, lx2, lx3 = [nx // bs for nx, bs in zip(nx_vals, block_size)]

    # Set coordinate system and related functions
    # coord = "cartesian"  # Adjust based on your data
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    # Define center functions if not provided
    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    # Populate coordinate arrays
    center_funcs = [center_func_1, center_func_2, center_func_3]
    for d in range(1, 4):
        xf = f"x{d}f"
        xv = f"x{d}v"
        nx = nx_vals[d - 1]
        if nx == 1:
            xmin = filedata[f"x{d}min"]
            xmax = filedata[f"x{d}max"]
            data[xf] = np.array([xmin, xmax], dtype=dtype)
        else:
            xmin = filedata[f"x{d}min"]
            xmax = filedata[f"x{d}max"]
            data[xf] = np.linspace(xmin, xmax, nx + 1, dtype=dtype)
        data[xv] = np.empty(nx, dtype=dtype)
        for i in range(nx):
            data[xv][i] = center_funcs[d - 1](data[xf][i], data[xf][i + 1])

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max = 0, nx1
    j_min, j_max = 0, nx2
    k_min, k_max = 0, nx3
    if x1_min is not None:
        i_min = max(i_min, np.searchsorted(data["x1f"], x1_min))
    if x1_max is not None:
        i_max = min(i_max, np.searchsorted(data["x1f"], x1_max))
    if x2_min is not None:
        j_min = max(j_min, np.searchsorted(data["x2f"], x2_min))
    if x2_max is not None:
        j_max = min(j_max, np.searchsorted(data["x2f"], x2_max))
    if x3_min is not None:
        k_min = max(k_min, np.searchsorted(data["x3f"], x3_min))
    if x3_max is not None:
        k_max = min(k_max, np.searchsorted(data["x3f"], x3_max))

    # Prepare arrays for data and bookkeeping
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.empty(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)
    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Step 3: Process each block
    for block_num in range(filedata["n_mbs"]):
        block_level = levels[block_num]
        block_location = logical_locations[block_num]

        # Implement logic for prolongation, restriction, subsampling
        if block_level <= level:
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
                block_data = filedata["mb_data"][q][block_num]
                if s > 1:
                    block_data = np.repeat(
                        np.repeat(np.repeat(block_data, s, axis=2), s, axis=1),
                        s,
                        axis=0,
                    )
                data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_data[
                    kl_s:ku_s, jl_s:ju_s, il_s:iu_s
                ]
        else:
            # Implement restriction logic here (similar to athdf function)
            pass

        if return_levels:
            data["Levels"][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level

    # Step 4: Finalize data
    if level < max_level and not subsample and not fast_restrict:
        # Remove volume factors from restricted data
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
                        for k in range(kl, ku):
                            for j in range(jl, ju):
                                for i in range(il, iu):
                                    x1m, x1p = data["x1f"][i], data["x1f"][i + 1]
                                    x2m, x2p = data["x2f"][j], data["x2f"][j + 1]
                                    x3m, x3p = data["x3f"][k], data["x3f"][k + 1]
                                    vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                    for q in quantities:
                                        data[q][k, j, i] /= vol

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = max_level

    return data


def read_single_rank_binary_as_athdf(
    filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    return_levels=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    *,
    meshblock_index_in_file=0,
):
    """
    Reads a single rank binary file and organizes data similar to
    athdf format without writing to file.
    """
    # Step 1: Read binary data for a single rank
    filedata = read_binary(filename)

    if isinstance(meshblock_index_in_file, (bool, np.bool_)) or not isinstance(
        meshblock_index_in_file, Integral
    ):
        raise TypeError("meshblock_index_in_file must be an integer")
    meshblock_index_in_file = int(meshblock_index_in_file)

    if raw:
        if meshblock_index_in_file != 0:
            raise ValueError("meshblock_index_in_file must be 0 when raw=True")
        return filedata

    _require_meshblocks_for_athdf(filedata, filename)
    if not 0 <= meshblock_index_in_file < filedata["n_mbs"]:
        raise IndexError(
            f"meshblock_index_in_file {meshblock_index_in_file} is out of range "
            f"for {filedata['n_mbs']} meshblocks"
        )

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    block_size = [filedata["nx1_mb"], filedata["nx2_mb"], filedata["nx3_mb"]]
    if dtype is None:
        dtype = np.float32

    # Set coordinate system and related functions
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

    # Populate coordinate arrays
    center_funcs = [center_func_1, center_func_2, center_func_3]
    for d in range(1, 4):
        xf = f"x{d}f"
        xv = f"x{d}v"
        nx = block_size[d - 1]

        # Use the meshblock geometry for local min and max
        xmin = filedata["mb_geometry"][meshblock_index_in_file, (d - 1) * 2]
        xmax = filedata["mb_geometry"][meshblock_index_in_file, (d - 1) * 2 + 1]

        data[xf] = np.linspace(xmin, xmax, nx + 1, dtype=dtype)
        data[xv] = np.empty(nx, dtype=dtype)
        for i in range(nx):
            data[xv][i] = center_funcs[d - 1](data[xf][i], data[xf][i + 1])

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max = 0, block_size[0]
    j_min, j_max = 0, block_size[1]
    k_min, k_max = 0, block_size[2]
    if x1_min is not None:
        i_min = max(i_min, np.searchsorted(data["x1f"], x1_min))
    if x1_max is not None:
        i_max = min(i_max, np.searchsorted(data["x1f"], x1_max))
    if x2_min is not None:
        j_min = max(j_min, np.searchsorted(data["x2f"], x2_min))
    if x2_max is not None:
        j_max = min(j_max, np.searchsorted(data["x2f"], x2_max))
    if x3_min is not None:
        k_min = max(k_min, np.searchsorted(data["x3f"], x3_min))
    if x3_max is not None:
        k_max = min(k_max, np.searchsorted(data["x3f"], x3_max))

    # Prepare arrays for data
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.empty(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)

    # Process the single block
    for q in quantities:
        block_data = filedata["mb_data"][q][meshblock_index_in_file]
        data[q] = block_data[k_min:k_max, j_min:j_max, i_min:i_max]

    if return_levels:
        data["Levels"].fill(filedata["mb_logical"][meshblock_index_in_file, 3])

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = filedata["mb_logical"][meshblock_index_in_file, 3]

    return data


def read_coarsened_binary_as_athdf(
    filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
):
    """
    Reads a bin file and organizes data similar to athdf format without writing to file.
    """
    # Step 1: Read binary data
    filedata = read_coarsened_binary(filename)

    # Step 2: Organize data similar to athdf
    if raw:
        return filedata

    _require_meshblocks_for_athdf(filedata, filename)

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    max_level = max(filedata["mb_logical"][:, 3])
    if level is None:
        level = max_level
    block_size = [filedata["nx1_mb"], filedata["nx2_mb"], filedata["nx3_mb"]]
    root_grid_size = [filedata["Nx1"], filedata["Nx2"], filedata["Nx3"]]
    levels = filedata["mb_logical"][:, 3]
    logical_locations = filedata["mb_logical"][:, :3]
    if dtype is None:
        dtype = np.float32

    # Calculate nx_vals
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:
            # Implement logic for sum or slice as in athdf
            nx_vals.append(root_grid_size[d] * 2**level)
        elif block_size[d] == 1:
            nx_vals.append(1)
        else:
            nx_vals.append(root_grid_size[d] * 2**level + 2 * num_ghost)
    nx1, nx2, nx3 = nx_vals
    lx1, lx2, lx3 = [nx // bs for nx, bs in zip(nx_vals, block_size)]

    # Set coordinate system and related functions
    # coord = "cartesian"  # Adjust based on your data
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    # Define center functions if not provided
    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    # Populate coordinate arrays
    center_funcs = [center_func_1, center_func_2, center_func_3]
    for d in range(1, 4):
        xf = f"x{d}f"
        xv = f"x{d}v"
        nx = nx_vals[d - 1]
        if nx == 1:
            xmin = filedata[f"x{d}min"]
            xmax = filedata[f"x{d}max"]
            data[xf] = np.array([xmin, xmax], dtype=dtype)
        else:
            xmin = filedata[f"x{d}min"]
            xmax = filedata[f"x{d}max"]
            data[xf] = np.linspace(xmin, xmax, nx + 1, dtype=dtype)
        data[xv] = np.empty(nx, dtype=dtype)
        for i in range(nx):
            data[xv][i] = center_funcs[d - 1](data[xf][i], data[xf][i + 1])

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max = 0, nx1
    j_min, j_max = 0, nx2
    k_min, k_max = 0, nx3
    if x1_min is not None:
        i_min = max(i_min, np.searchsorted(data["x1f"], x1_min))
    if x1_max is not None:
        i_max = min(i_max, np.searchsorted(data["x1f"], x1_max))
    if x2_min is not None:
        j_min = max(j_min, np.searchsorted(data["x2f"], x2_min))
    if x2_max is not None:
        j_max = min(j_max, np.searchsorted(data["x2f"], x2_max))
    if x3_min is not None:
        k_min = max(k_min, np.searchsorted(data["x3f"], x3_min))
    if x3_max is not None:
        k_max = min(k_max, np.searchsorted(data["x3f"], x3_max))

    # Prepare arrays for data and bookkeeping
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.empty(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)
    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Step 3: Process each block
    for block_num in range(filedata["n_mbs"]):
        block_level = levels[block_num]
        block_location = logical_locations[block_num]

        # Implement logic for prolongation, restriction, subsampling
        if block_level <= level:
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
                block_data = filedata["mb_data"][q][block_num]
                if s > 1:
                    block_data = np.repeat(
                        np.repeat(np.repeat(block_data, s, axis=2), s, axis=1),
                        s,
                        axis=0,
                    )
                data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_data[
                    kl_s:ku_s, jl_s:ju_s, il_s:iu_s
                ]
        else:
            # Implement restriction logic here (similar to athdf function)
            pass

        if return_levels:
            data["Levels"][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level

    # Step 4: Finalize data
    if level < max_level and not subsample and not fast_restrict:
        # Remove volume factors from restricted data
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
                        for k in range(kl, ku):
                            for j in range(jl, ju):
                                for i in range(il, iu):
                                    x1m, x1p = data["x1f"][i], data["x1f"][i + 1]
                                    x2m, x2p = data["x2f"][j], data["x2f"][j + 1]
                                    x3m, x3p = data["x3f"][k], data["x3f"][k + 1]
                                    vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                    for q in quantities:
                                        data[q][k, j, i] /= vol

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = max_level

    return data


def write_athdf(filename, fdata, varsize_bytes=4, locsize_bytes=8):
    """
    Writes an athdf (hdf5) file from a loaded python filedata object.

    args:
      filename      - string
          filename for output athdf (hdf5) file
      fdata         - dict
          dictionary of fluid file data, e.g., as loaded from read_binary(...)
      varsize_bytes - int (default=4, options=4,8)
          number of bytes to use for output variable data
      locsize_bytes - int (default=8, options=4,8)
          number of bytes to use for output location data
    """

    if varsize_bytes not in [4, 8]:
        raise ValueError(f"varsizebytes must be 4 or 8, not {varsize_bytes}")
    if locsize_bytes not in [4, 8]:
        raise ValueError(f"locsizebytes must be 4 or 8, not {locsize_bytes}")
    locfmt = "<f4" if locsize_bytes == 4 else "<f8"
    varfmt = "<f4" if varsize_bytes == 4 else "<f8"

    # extract Mesh/MeshBlock parameters
    nmb = fdata["n_mbs"]
    Nx1 = fdata["Nx1"]  # noqa: F841
    Nx2 = fdata["Nx2"]
    Nx3 = fdata["Nx3"]
    nx1 = fdata["nx1_mb"]
    nx2 = fdata["nx2_mb"]
    nx3 = fdata["nx3_mb"]
    nx1_out = fdata["nx1_out_mb"]
    nx2_out = fdata["nx2_out_mb"]
    nx3_out = fdata["nx3_out_mb"]

    number_of_moments = fdata.get("number_of_moments", 1)

    # check dimensionality/slicing
    two_d = Nx2 != 1 and Nx3 == 1
    three_d = Nx3 != 1
    x1slice = nx1_out == 1
    x2slice = nx2_out == 1 and (two_d or three_d)
    x3slice = nx3_out == 1 and three_d

    # keep variable order but separate out magnetic field
    vars_without_b = [v for v in fdata["var_names"] if "bcc" not in v]
    vars_only_b = [v for v in fdata["var_names"] if v not in vars_without_b]

    if len(vars_only_b) > 0:
        B = np.zeros((3 * number_of_moments, nmb, nx3_out, nx2_out, nx1_out))
    Levels = np.zeros(nmb)
    LogicalLocations = np.zeros((nmb, 3))
    uov = np.zeros((len(vars_without_b), nmb, nx3_out, nx2_out, nx1_out))
    x1f = np.zeros((nmb, nx1_out + 1))
    x1v = np.zeros((nmb, nx1_out))
    x2f = np.zeros((nmb, nx2_out + 1))
    x2v = np.zeros((nmb, nx2_out))
    x3f = np.zeros((nmb, nx3_out + 1))
    x3v = np.zeros((nmb, nx3_out))

    for ivar, var in enumerate(vars_without_b):
        uov[ivar] = fdata["mb_data"][var]
    for ibvar, bvar in enumerate(vars_only_b):
        B[ibvar] = fdata["mb_data"][bvar]

    for mb in range(nmb):
        logical = fdata["mb_logical"][mb]
        LogicalLocations[mb] = logical[:3]
        Levels[mb] = logical[-1]
        geometry = fdata["mb_geometry"][mb]
        mb_x1f = np.linspace(geometry[0], geometry[1], nx1 + 1)
        mb_x1v = 0.5 * (mb_x1f[1:] + mb_x1f[:-1])
        mb_x2f = np.linspace(geometry[2], geometry[3], nx2 + 1)
        mb_x2v = 0.5 * (mb_x2f[1:] + mb_x2f[:-1])
        mb_x3f = np.linspace(geometry[4], geometry[5], nx3 + 1)
        mb_x3v = 0.5 * (mb_x3f[1:] + mb_x3f[:-1])
        if x1slice:
            x1f[mb] = np.array(
                mb_x1f[(fdata["mb_index"][mb][0]): (fdata["mb_index"][mb][0] + 2)]
            )
            x1v[mb] = np.array([np.average(mb_x1f)])
        else:
            x1f[mb] = mb_x1f
            x1v[mb] = mb_x1v
        if x2slice:
            x2f[mb] = np.array(
                mb_x2f[(fdata["mb_index"][mb][2]): (fdata["mb_index"][mb][2] + 2)]
            )
            x2v[mb] = np.array([np.average(x2f[mb])])
        else:
            x2f[mb] = mb_x2f
            x2v[mb] = mb_x2v
        if x3slice:
            x3f[mb] = np.array(
                mb_x3f[(fdata["mb_index"][mb][4]): (fdata["mb_index"][mb][4] + 2)]
            )
            x3v[mb] = np.array([np.average(x3f[mb])])
        else:
            x3f[mb] = mb_x3f
            x3v[mb] = mb_x3v

    # set dataset names and number of variables
    dataset_names = [np.array("uov", dtype="|S21")]
    dataset_nvars = [len(vars_without_b)]
    if len(vars_only_b) > 0:
        dataset_names.append(np.array("B", dtype="|S21"))
        dataset_nvars.append(len(vars_only_b))

    # Set Attributes
    hfp = h5py.File(filename, "w")
    hfp.attrs["Header"] = fdata["header"]
    hfp.attrs["Time"] = fdata["time"]
    hfp.attrs["NumCycles"] = fdata["cycle"]
    hfp.attrs["Coordinates"] = np.array("cartesian", dtype="|S11")
    hfp.attrs["NumMeshBlocks"] = fdata["n_mbs"]
    hfp.attrs["MaxLevel"] = int(max(Levels))
    hfp.attrs["MeshBlockSize"] = [
        fdata["nx1_out_mb"],
        fdata["nx2_out_mb"],
        fdata["nx3_out_mb"],
    ]
    hfp.attrs["RootGridSize"] = [fdata["Nx1"], fdata["Nx2"], fdata["Nx3"]]
    hfp.attrs["RootGridX1"] = [fdata["x1min"], fdata["x1max"], 1.0]
    hfp.attrs["RootGridX2"] = [fdata["x2min"], fdata["x2max"], 1.0]
    hfp.attrs["RootGridX3"] = [fdata["x3min"], fdata["x3max"], 1.0]
    hfp.attrs["DatasetNames"] = dataset_names
    hfp.attrs["NumVariables"] = dataset_nvars
    hfp.attrs["VariableNames"] = [
        np.array(i, dtype="|S21") for i in (vars_without_b + vars_only_b)
    ]

    # Create Datasets
    if len(vars_only_b) > 0:
        hfp.create_dataset("B", data=B, dtype=varfmt)
    hfp.create_dataset("Levels", data=Levels, dtype=">i4")
    hfp.create_dataset("LogicalLocations", data=LogicalLocations, dtype=">i8")
    hfp.create_dataset("uov", data=uov, dtype=varfmt)
    hfp.create_dataset("x1f", data=x1f, dtype=locfmt)
    hfp.create_dataset("x1v", data=x1v, dtype=locfmt)
    hfp.create_dataset("x2f", data=x2f, dtype=locfmt)
    hfp.create_dataset("x2v", data=x2v, dtype=locfmt)
    hfp.create_dataset("x3f", data=x3f, dtype=locfmt)
    hfp.create_dataset("x3v", data=x3v, dtype=locfmt)
    hfp.close()


def write_xdmf_for(xdmfname, dumpname, fdata, mode="auto"):
    """
    Writes an xdmf file for a fluid snapshot file.

    args:
      xdmfname - string
          name of xdmf file
      dumpname - string
          location of fluid data file relative to xdmfname directory
      fdata    - dict
          dictionary of fluid file data, e.g., as loaded from read_binary(...)
      mode     - string (unimplemented)
          force xdmf for format (auto sets by extension)
    """

    fp = open(xdmfname, "w")

    def write_meshblock(fp, mb, nx1, nx2, nx3, nmb, dumpname, vars_no_b, vars_w_b):
        fp.write(f"""  <Grid Name="MeshBlock{mb}" GridType="Uniform">\n""")
        fp.write("""   <Topology TopologyType="3DRectMesh" """)
        fp.write(f""" NumberOfElements="{nx3+1} {nx2+1} {nx1+1}"/>\n""")
        fp.write("""   <Geometry GeometryType="VXVYVZ">\n""")
        fp.write(
            f"""    <DataItem ItemType="HyperSlab" Dimensions="{nx1+1}">
     <DataItem Dimensions="3 2" NumberType="Int"> {mb} 0 1 1 1 {nx1+1} </DataItem>
     <DataItem Dimensions="{nmb} {nx1+1}" Format="HDF"> {dumpname}:/x1f </DataItem>
    </DataItem>
    <DataItem ItemType="HyperSlab" Dimensions="{nx2+1}">
     <DataItem Dimensions="3 2" NumberType="Int"> {mb} 0 1 1 1 {nx2+1} </DataItem>
     <DataItem Dimensions="{nmb} {nx2+1}" Format="HDF"> {dumpname}:/x2f </DataItem>
    </DataItem>
    <DataItem ItemType="HyperSlab" Dimensions="{nx3+1}">
     <DataItem Dimensions="3 2" NumberType="Int"> {mb} 0 1 1 1 {nx3+1} </DataItem>
     <DataItem Dimensions="{nmb} {nx3+1}" Format="HDF"> {dumpname}:/x3f </DataItem>
    </DataItem>
   </Geometry>\n"""
        )

        nvar_no_b = len(vars_no_b)
        for vi, var_name in enumerate(vars_no_b):
            fp.write(
                f"""   <Attribute Name="{var_name}" Center="Cell">
    <DataItem ItemType="HyperSlab" Dimensions="{nx3} {nx2} {nx1}">
     <DataItem Dimensions="3 5" NumberType="Int">
      {vi} {mb} 0 0 0 1 1 1 1 1 1 1 {nx3} {nx2} {nx1}
     </DataItem>
     <DataItem Dimensions="{nvar_no_b} {nmb} {nx3} {nx2} {nx1}" Format="HDF">
      {dumpname}:/uov
     </DataItem>
    </DataItem>
   </Attribute>\n"""
            )

        nvar_w_b = len(vars_w_b)
        if nvar_w_b > 0:
            for vi, var_name in enumerate(vars_w_b):
                fp.write(
                    f"""   <Attribute Name="{var_name}" Center="Cell">
        <DataItem ItemType="HyperSlab" Dimensions="{nx3} {nx2} {nx1}">
         <DataItem Dimensions="3 5" NumberType="Int">
          {vi} {mb} 0 0 0 1 1 1 1 1 1 1 {nx3} {nx2} {nx1}
         </DataItem>
         <DataItem Dimensions="{nvar_w_b} {nmb} {nx3} {nx2} {nx1}" Format="HDF">
          {dumpname}:/B
         </DataItem>
        </DataItem>
       </Attribute>\n"""
                )

        fp.write("""  </Grid>\n""")

    fp.write(
        """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
<Information Name="TimeVaryingMetaData" Value="True"/>\n"""
    )
    fp.write("""<Domain>\n""")
    fp.write("""<Grid Name="Mesh" GridType="Collection">\n""")
    fp.write(f""" <Time Value="{fdata['time']}"/>\n""")

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

    fp.write("""</Grid>\n""")
    fp.write("""</Domain>\n""")
    fp.write("""</Xdmf>\n""")

    fp.close()


def convert_file(binary_fname, assemble_shards=False, coarsened=None):
    """
    Converts a binary file, optionally assembling its rank/node shard family.

    args:
      binary_filename - string
        filename of bin file to convert
      assemble_shards - bool, optional
        assemble sibling rank_* or node_* files before writing output
      coarsened - bool or None, optional
        use the coarsened-binary reader; when None, infer from a .cbin suffix

    This will create new files "binary_data.bin" -> "binary_data.athdf" and
    "binary_data.athdf.xdmf"
    """
    athdf_fname = os.path.splitext(binary_fname)[0] + ".athdf"
    xdmf_fname = athdf_fname + ".xdmf"
    if coarsened is None:
        coarsened = binary_fname.endswith(".cbin")
    if coarsened:
        filedata = read_coarsened_binary(binary_fname, assemble_shards=assemble_shards)
    else:
        filedata = read_binary(binary_fname, assemble_shards=assemble_shards)
    if filedata["n_mbs"] == 0:
        raise ValueError(
            f"cannot convert {binary_fname!r}: binary output contains no meshblocks"
        )
    write_athdf(athdf_fname, filedata)
    write_xdmf_for(xdmf_fname, os.path.basename(athdf_fname), filedata)


__all__ = [
    "read_binary",
    "read_coarsened_binary",
    "read_all_ranks_binary",
    "read_all_ranks_coarsened_binary",
    "read_binary_as_athdf",
    "read_all_ranks_binary_as_athdf",
    "read_all_ranks_coarsened_binary_as_athdf",
    "read_single_rank_binary_as_athdf",
    "read_rank_binary_as_athdf",
    "read_coarsened_binary_as_athdf",
    "write_athdf",
    "write_xdmf_for",
    "convert_file",
]


if __name__ == "__main__":
    import argparse
    import sys

    try:
        from tqdm import tqdm
    except ModuleNotFoundError:

        def tqdm(L):
            for x in L:
                print(x)
                yield x

    parser = argparse.ArgumentParser(
        description="Convert AthenaK binary output to ATHDF/XDMF."
    )
    parser.add_argument(
        "--assemble-shards",
        action="store_true",
        help="assemble matching rank_* or node_* sibling shards before conversion",
    )
    parser.add_argument(
        "--coarsened",
        action="store_true",
        default=None,
        help="force the coarsened-binary reader (otherwise inferred from .cbin)",
    )
    parser.add_argument("binary_files", nargs="+", help="binary files to convert")
    args = parser.parse_args(sys.argv[1:])

    for binary_fname in tqdm(args.binary_files):
        convert_file(
            binary_fname,
            assemble_shards=args.assemble_shards,
            coarsened=args.coarsened,
        )
