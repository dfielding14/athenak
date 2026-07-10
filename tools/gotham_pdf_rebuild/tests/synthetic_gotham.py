#!/usr/bin/env python3
"""Generate small, deterministic AthenaK v1.1 GOTHAM hydro_w shards."""

from __future__ import annotations

import argparse
import json
import struct
import sys
from array import array
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterator, Sequence, Tuple


DEFAULT_NX = 32
VARIABLES = ("dens", "velx", "vely", "velz", "eint", "s_00", "s_01", "s_02")
SEQUENCE_DEFAULT = "00000"
TIME_DEFAULT = 12.5
CYCLE_DEFAULT = 42
GAMMA_DEFAULT = 5.0 / 3.0

_RHO = (1.0e-8, 1.0e-5, 1.0e-2, 1.0, 1.0e3, 1.0e5, 1.0e8, 0.25)
_TEMPERATURE = (1.0e-6, 1.0e-3, 1.0e-1, 1.0, 1.0e2, 1.0e4, 7.0e4, 1.0e6)
_VX = (-1000.0, -100.0, -10.0, -1.0, 0.0, 1.0, 100.0, 1000.0)
_VY = (500.0, -50.0, 5.0, -0.5, 0.0, 0.5, 50.0, -500.0)
_VZ = (-200.0, 20.0, -2.0, 0.0, 0.0, 2.0, -20.0, 200.0)
_S0 = (1.0e-8, 1.0e-5, 1.0e-4, 1.0e-2, 0.1, 0.3, 1.0, 0.05)
_S1 = (1.0e-8, 1.0e-5, 1.0e-4, 1.0e-2, 0.1, 1.0, 2.0, 0.02)
_S2 = (2.0, 1.0, 0.1, 0.01, 1.0e-4, 1.0e-5, 1.0e-8, 0.03)


@dataclass(frozen=True)
class BlockSpec:
    """One deterministic synthetic MeshBlock."""

    global_id: int
    geometry: Tuple[float, float, float, float, float, float]
    shape: Tuple[int, int, int]

    @property
    def cells(self) -> int:
        return self.shape[0] * self.shape[1] * self.shape[2]


def f32(value: float) -> float:
    """Round exactly as AthenaK's float binary output does."""

    return struct.unpack("=f", struct.pack("=f", value))[0]


def make_blocks(count: int, shape: Tuple[int, int, int]) -> list[BlockSpec]:
    if count <= 0:
        raise ValueError("block count must be positive")
    if any(extent <= 0 for extent in shape):
        raise ValueError("MeshBlock extents must be positive")
    domain_min = -256.0
    domain_max = 256.0

    def left_edge(index: int) -> float:
        fraction = index / count
        return (
            (fraction * domain_max - fraction * domain_min)
            - (0.5 * domain_max - 0.5 * domain_min)
            + (0.5 * domain_min + 0.5 * domain_max)
        )

    return [
        BlockSpec(
            global_id=index,
            geometry=(
                domain_min if index == 0 else left_edge(index),
                domain_max if index == count - 1 else left_edge(index + 1),
                domain_min,
                domain_max,
                domain_min,
                domain_max,
            ),
            shape=shape,
        )
        for index in range(count)
    ]


def field_values(block_id: int, cell: int) -> Tuple[float, ...]:
    """Return the eight stored float32 values for one cell."""

    rho = f32(_RHO[(cell + block_id) % len(_RHO)])
    temperature = _TEMPERATURE[(cell // 3 + 2 * block_id) % len(_TEMPERATURE)]
    vx = _VX[(cell // 5 + block_id) % len(_VX)]
    vy = _VY[(cell // 7 + 3 * block_id) % len(_VY)]
    vz = _VZ[(cell // 11 + 5 * block_id) % len(_VZ)]
    return (
        rho,
        f32(vx),
        f32(vy),
        f32(vz),
        f32(rho * temperature),
        f32(_S0[(cell // 13 + block_id) % len(_S0)]),
        f32(_S1[(cell // 17 + 2 * block_id) % len(_S1)]),
        f32(_S2[(cell // 19 + 3 * block_id) % len(_S2)]),
    )


def cell_center(block: BlockSpec, cell: int) -> Tuple[float, float, float, float]:
    nx1, nx2, nx3 = block.shape
    i = cell % nx1
    j = (cell // nx1) % nx2
    k = cell // (nx1 * nx2)
    x0, x1, y0, y1, z0, z1 = block.geometry
    dx = (x1 - x0) / nx1
    dy = (y1 - y0) / nx2
    dz = (z1 - z0) / nx3
    return x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy, z0 + (k + 0.5) * dz, dx * dy * dz


def iter_cells(blocks: Sequence[BlockSpec]) -> Iterator[Tuple[BlockSpec, int, Tuple[float, ...]]]:
    for block in blocks:
        for cell in range(block.cells):
            yield block, cell, field_values(block.global_id, cell)


def _embedded_input(
    gamma: float,
    shape: Tuple[int, int, int],
    root_blocks: Tuple[int, int, int],
    domain_min: Tuple[float, float, float],
    domain_max: Tuple[float, float, float],
) -> bytes:
    nx1, nx2, nx3 = shape
    nmb1, nmb2, nmb3 = root_blocks
    return (
        "<mesh>\n"
        f"nx1 = {nx1 * nmb1}\n"
        f"x1min = {domain_min[0]:.17g}\n"
        f"x1max = {domain_max[0]:.17g}\n"
        f"nx2 = {nx2 * nmb2}\n"
        f"x2min = {domain_min[1]:.17g}\n"
        f"x2max = {domain_max[1]:.17g}\n"
        f"nx3 = {nx3 * nmb3}\n"
        f"x3min = {domain_min[2]:.17g}\n"
        f"x3max = {domain_max[2]:.17g}\n"
        "\n"
        "<meshblock>\n"
        f"nx1 = {nx1}\n"
        f"nx2 = {nx2}\n"
        f"nx3 = {nx3}\n"
        "\n<hydro>\n"
        f"gamma = {gamma:.17g}\n"
        "cgm_cooling = true\n"
        "hrate = 2e-26\n"
        "hscale_norm = 1347.6541799847946\n"
        "hscale_radius = 0.9402361259893334\n"
        "hscale_height = 0.04950000000848683\n"
        "T_max = 5e8\n"
        "\n<units>\n"
        "length_cgs = 3.0856775809623245e+21\n"
        "mass_cgs = 3.036951775493658e+39\n"
        "time_cgs = 3.15576e+15\n"
        "mu = 0.62\n"
    ).encode("ascii")


def _header(
    time: float,
    cycle: int,
    gamma: float,
    shape: Tuple[int, int, int],
    root_blocks: Tuple[int, int, int],
    domain_min: Tuple[float, float, float],
    domain_max: Tuple[float, float, float],
) -> bytes:
    embedded = _embedded_input(
        gamma, shape, root_blocks, domain_min, domain_max
    )
    prefix = (
        "Athena binary output version=1.1\n"
        "  size of preheader=5\n"
        f"  time={time:.17g}\n"
        f"  cycle={cycle}\n"
        "  size of location=8\n"
        "  size of variable=4\n"
        f"  number of variables={len(VARIABLES)}\n"
        f"  variables:  {'  '.join(VARIABLES)}  \n"
        f"  header offset={len(embedded)}\n"
    ).encode("ascii")
    return prefix + embedded


def _write_block(stream, block: BlockSpec) -> None:
    nx1, nx2, nx3 = block.shape
    ints = (0, nx1 - 1, 0, nx2 - 1, 0, nx3 - 1, block.global_id, 0, 0, 0)
    stream.write(struct.pack("=10i6d", *(ints + block.geometry)))

    fields = [array("f") for _ in VARIABLES]
    for cell in range(block.cells):
        for destination, value in zip(fields, field_values(block.global_id, cell)):
            destination.append(value)
    for values in fields:
        stream.write(values.tobytes())


def write_dataset(
    output_dir: Path,
    *,
    blocks: int = 4,
    shards: int = 2,
    sequence: str = SEQUENCE_DEFAULT,
    time: float = TIME_DEFAULT,
    cycle: int = CYCLE_DEFAULT,
    gamma: float = GAMMA_DEFAULT,
    nx1: int = DEFAULT_NX,
    nx2: int | None = None,
    nx3: int | None = None,
) -> list[BlockSpec]:
    """Write a synthetic per-node shard dataset and return its block contract."""

    if sys.byteorder != "little":
        raise RuntimeError("The synthetic AthenaK writer currently requires little endian")
    if shards <= 0:
        raise ValueError("shard count must be positive")
    if len(sequence) != 5 or not sequence.isdigit():
        raise ValueError("sequence must contain exactly five digits")

    shape = (nx1, nx1 if nx2 is None else nx2, nx1 if nx3 is None else nx3)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    block_specs = make_blocks(blocks, shape)
    assignments = [[] for _ in range(shards)]
    for index, block in enumerate(block_specs):
        assignments[index % shards].append(block)

    domain_min = (
        min(block.geometry[0] for block in block_specs),
        min(block.geometry[2] for block in block_specs),
        min(block.geometry[4] for block in block_specs),
    )
    domain_max = (
        max(block.geometry[1] for block in block_specs),
        max(block.geometry[3] for block in block_specs),
        max(block.geometry[5] for block in block_specs),
    )
    header = _header(
        time,
        cycle,
        gamma,
        shape,
        (blocks, 1, 1),
        domain_min,
        domain_max,
    )
    for shard, assigned in enumerate(assignments):
        shard_dir = output_dir / f"node_{shard:08d}"
        shard_dir.mkdir(parents=True, exist_ok=True)
        path = shard_dir / f"gotham.hydro_w.{sequence}.bin"
        with path.open("wb") as stream:
            stream.write(header)
            for block in assigned:
                _write_block(stream, block)

    manifest = {
        "format": "Athena binary output version=1.1",
        "sequence": sequence,
        "time": time,
        "cycle": cycle,
        "gamma": gamma,
        "meshblock_shape": list(shape),
        "variables": list(VARIABLES),
        "blocks": [asdict(block) for block in block_specs],
        "shards": shards,
    }
    (output_dir / "synthetic_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return block_specs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--blocks", type=int, default=4)
    parser.add_argument("--shards", type=int, default=2)
    parser.add_argument("--sequence", default=SEQUENCE_DEFAULT)
    parser.add_argument("--time", type=float, default=TIME_DEFAULT)
    parser.add_argument("--cycle", type=int, default=CYCLE_DEFAULT)
    parser.add_argument("--gamma", type=float, default=GAMMA_DEFAULT)
    parser.add_argument("--nx1", type=int, default=DEFAULT_NX)
    parser.add_argument("--nx2", type=int)
    parser.add_argument("--nx3", type=int)
    args = parser.parse_args()
    blocks = write_dataset(
        args.output_dir,
        blocks=args.blocks,
        shards=args.shards,
        sequence=args.sequence,
        time=args.time,
        cycle=args.cycle,
        gamma=args.gamma,
        nx1=args.nx1,
        nx2=args.nx2,
        nx3=args.nx3,
    )
    payload_bytes = sum(10 * 4 + 6 * 8 + len(VARIABLES) * block.cells * 4 for block in blocks)
    print(
        f"wrote {len(blocks)} blocks across {args.shards} shards "
        f"({payload_bytes / 1024**2:.3f} MiB records)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
