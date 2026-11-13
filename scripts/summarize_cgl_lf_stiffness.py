"""Summarize one cycle/phase of rank-local LF maps using one meshblock at a time.

Usage: python3 scripts/summarize_cgl_lf_stiffness.py lf-stiffness.123.pre.rank*.bin
Nearest-rank median/p99 are bracketed by bins 1/32 octave wide (about 2.2%).
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np


def read_header(path):
    path = Path(path)
    with path.open("rb") as stream:
        if stream.readline() != b"# AthenaK LF stiffness v1\n":
            raise ValueError(f"{path}: unsupported LF map header")
        pairs = [word.split("=", 1) for word in stream.readline().decode().split()]
        if any(len(pair) != 2 for pair in pairs) or len(dict(pairs)) != len(pairs):
            raise ValueError(f"{path}: invalid or duplicate metadata")
        meta = dict(pairs)
        for key in ("cycle", "rank", "nranks", "nblocks", "total_blocks", "nx1",
                    "nx2", "nx3", "real_bytes", "nstages"):
            meta[key] = int(meta[key])
        if meta["real_bytes"] not in (4, 8):
            raise ValueError(f"{path}: unsupported LF map precision")
        real = np.dtype(f"f{meta['real_bytes']}").type
        for key in ("time", "dt_outer", "dt_sweep", "sts_safety", "dt_sts",
                    "raw_dt_rank"):
            meta[key] = float(real(meta[key]))
        if (meta["endian"] not in ("little", "big")
                or meta["phase"] not in ("pre", "post")
                or not 0 <= meta["rank"] < meta["nranks"]
                or meta["cycle"] < 0 or meta["nstages"] < 0
                or any(meta[key] <= 0 for key in
                       ("nblocks", "total_blocks", "nx1", "nx2", "nx3"))
                or any(not math.isfinite(meta[key]) or meta[key] <= 0 for key in
                       ("dt_outer", "dt_sweep", "sts_safety", "dt_sts"))
                or meta["sts_safety"] > 1 or not math.isfinite(meta["time"])
                or not math.isfinite(meta["raw_dt_rank"]) or meta["raw_dt_rank"] < 0):
            raise ValueError(f"{path}: invalid LF map metadata")
        blocks = []
        for _ in range(meta["nblocks"]):
            row = stream.readline().split()
            if len(row) != 6 or row[0] != b"block":
                raise ValueError(f"{path}: invalid block header")
            blocks.append(tuple(map(int, row[1:])))
        if stream.readline() != b"data\n":
            raise ValueError(f"{path}: missing data marker")
        offset = stream.tell()
    shape = (meta["nblocks"], meta["nx3"], meta["nx2"], meta["nx1"])
    count = math.prod(shape)
    if path.stat().st_size != offset + count * meta["real_bytes"]:
        raise ValueError(f"{path}: truncated or oversized LF map payload")
    return meta, blocks, offset


def rkl_stages(ratio):
    """Smallest sufficient odd stage count for an active RKL2 sweep."""
    if ratio == 0:
        return 0
    maximum = 2**31 - 1
    if not math.isfinite(ratio) or ratio > (maximum*maximum + maximum - 2)/4:
        return "exceeds int32"
    stages = max(3, math.ceil((-1 + math.sqrt(9 + 16*ratio))/2))
    stages += stages % 2 == 0
    while (stages*stages + stages - 2)/4 < ratio:
        stages += 2
    while stages > 3 and ((stages-2)**2 + stages-2 - 2)/4 >= ratio:
        stages -= 2
    return stages


def stiffness(log_value):
    if log_value == -math.inf:
        return 0.0
    if log_value >= 1024:
        return "infinite" if log_value == math.inf else f"2^{log_value:.8g}"
    return 2.0**log_value


def summarize(paths):
    if not paths:
        raise ValueError("provide all rank maps for one cycle and phase")
    maps = [(Path(path), *read_header(path)) for path in paths]
    first = maps[0][1]
    shared = set(first) - {"rank", "nblocks", "raw_dt_rank"}
    if any(set(meta) != set(first) or any(meta[key] != first[key] for key in shared)
           for _, meta, _, _ in maps):
        raise ValueError("mixed cycle, phase, timestep, geometry, or precision")
    ranks = [meta["rank"] for _, meta, _, _ in maps]
    if sorted(ranks) != list(range(first["nranks"])):
        raise ValueError("missing or duplicate rank maps")
    block_ids = [block[0] for _, _, blocks, _ in maps for block in blocks]
    if sorted(block_ids) != list(range(first["total_blocks"])):
        raise ValueError("missing or duplicate global meshblock IDs")

    log_scale = math.log2(first["dt_sweep"]) - math.log2(first["sts_safety"])
    histogram = np.zeros(2*4096*32, dtype=np.int64)
    total = inactive = invalid = 0
    cutoff = 10*min(meta["raw_dt_rank"] for _, meta, _, _ in maps)
    near_cells = near_blocks = 0
    locations = []
    for path, meta, blocks, offset in maps:
        dtype = np.dtype(("<" if meta["endian"] == "little" else ">") +
                         f"f{meta['real_bytes']}")
        shape = (meta["nx3"], meta["nx2"], meta["nx1"])
        rank_min = np.finfo(dtype).max
        with path.open("rb") as stream:
            stream.seek(offset)
            for block in blocks:
                data = np.fromfile(stream, dtype=dtype, count=math.prod(shape)).reshape(shape)
                if np.any(np.isnan(data) | (data < 0)):
                    raise ValueError(f"rank {meta['rank']}: NaN or negative local timestep")
                index = np.unravel_index(np.argmin(data), data.shape)
                minimum = float(data[index])
                rank_min = min(rank_min, minimum)
                invalid += int(np.count_nonzero(data == 0))
                inactive += int(np.count_nonzero(np.isposinf(data)))
                total += data.size
                count = np.count_nonzero(np.isfinite(data) & (data <= cutoff))
                near_cells += int(count)
                near_blocks += count > 0
                positive = data[(data > 0) & np.isfinite(data)].astype(np.float64)
                bins = np.floor((log_scale - np.log2(positive))*32).astype(np.int64)
                histogram += np.bincount(bins + 4096*32, minlength=histogram.size)
                locations.append(dict(rank=meta["rank"], gid=block[0], level=block[1],
                                      logical_location=list(block[2:]),
                                      cell_ijk=list(map(int, index[::-1])), raw_dt=minimum))
        if rank_min != meta["raw_dt_rank"]:
            raise ValueError(f"rank {meta['rank']}: map minimum differs from raw_dt_rank")

    locations.sort(key=lambda row: (row["raw_dt"], row["gid"]))
    worst_dt = locations[0]["raw_dt"]

    def log_ratio(dt):
        return math.inf if dt == 0 else log_scale - math.log2(dt)

    def percentile(fraction):
        position = max(1, math.ceil(fraction*total))
        if position <= inactive:
            return [0.0, 0.0]
        if position > total-invalid:
            return ["infinite", "infinite"]
        index = int(np.searchsorted(histogram.cumsum(), position-inactive))
        lower = index/32 - 4096
        return [stiffness(lower), stiffness(lower+1/32)]

    maximum = log_ratio(worst_dt)
    real = np.dtype(f"f{first['real_bytes']}").type

    def stage_count(raw_dt):
        # The controller multiplies in Real before promoting the ratio to long double.
        budget = real(first["sts_safety"])*real(raw_dt)
        return rkl_stages(np.longdouble(first["dt_sweep"])/budget if budget > 0 else math.inf)

    effective_dt = real(first["sts_safety"])*real(worst_dt)
    predicted = stage_count(worst_dt)
    next_stages = stage_count(locations[1]["raw_dt"]) if len(locations) > 1 else None
    warnings = []
    precision = np.finfo(f"f{first['real_bytes']}").eps
    if math.isfinite(effective_dt) and not math.isclose(
            effective_dt, first["dt_sts"], rel_tol=8*precision):
        warnings.append("selected dt_sts differs from sts_safety*minimum raw timestep")
    if predicted != first["nstages"]:
        warnings.append("fresh map predicts different stages from the selected sweep; "
                        "compare dt_sts with sts_safety*minimum raw timestep")
    # JSON uses strings for unbounded values, preserving valid strict JSON.
    for row in locations:
        if math.isinf(row["raw_dt"]):
            row["raw_dt"] = "inactive"
    return dict(cycle=first["cycle"], phase=first["phase"], time=first["time"],
                owned_cells=total, blocks=len(locations), invalid_bound_cells=invalid,
                inactive_cells=inactive, dt_outer=first["dt_outer"],
                dt_sweep=first["dt_sweep"], sts_safety=first["sts_safety"],
                dt_sts=first["dt_sts"], maximum_stiffness=stiffness(maximum),
                median_stiffness_interval=percentile(.5),
                p99_stiffness_interval=percentile(.99),
                cells_within_factor_ten=near_cells, blocks_within_factor_ten=near_blocks,
                worst_block=locations[0],
                next_worst_block=locations[1] if len(locations) > 1 else None,
                predicted_rkl_stages=predicted, selected_rkl_stages=first["nstages"],
                predicted_rkl_stages_without_worst_block=next_stages,
                warnings=warnings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("maps", nargs="+", type=Path)
    args = parser.parse_args()
    try:
        result = summarize(args.maps)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"LF stiffness summary failed: {error}\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
