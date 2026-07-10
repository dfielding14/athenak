#!/usr/bin/env python3
"""Add exact archived PDF output identities to the five static cube manifests."""

from __future__ import annotations

import csv
import math
from pathlib import Path


JOBS_DIR = Path(__file__).resolve().parent
MANIFEST_DIR = JOBS_DIR / "manifests"
CATALOG = Path(
    "/lustre/orion/ast207/proj-shared/gotham/analysis/catalog/generated/"
    "output_time_lookup.csv"
)
CONTROL_STREAM = "pdf_coord_costheta_coord_r_temperature_vel_sph_r"
ARCHIVE_ROOT = Path("/lustre/orion/ast207/proj-shared/brent/gotham")
OLD_FIELDS = (
    "simulation",
    "phase",
    "sequence",
    "input_relative",
    "expected_shards",
    "expected_nodes",
)
NEW_FIELDS = OLD_FIELDS + ("output_number", "output_time", "control_pdf_sequence")


def time_tolerance(value: float) -> float:
    if value == 0.0:
        return 0.5e-6
    exponent = math.floor(math.log10(abs(value)))
    return 0.51 * 10.0 ** (exponent - 5)


def read_catalog():
    cubes = {}
    controls = {}
    with CATALOG.open(newline="") as stream:
        for row in csv.DictReader(stream):
            key = (row["run"], row["phase"])
            if row["output_family"] == "binary" and row["stream"] == "hydro_w":
                cubes[(key, f"{int(row['sequence_index']):05d}")] = float(
                    row["time_code"]
                )
            elif row["output_family"] == "pdf" and row["stream"] == CONTROL_STREAM:
                controls.setdefault(key, []).append(
                    (f"{int(row['sequence_index']):05d}", float(row["time_code"]))
                )
    return cubes, controls


def read_cube_time(row) -> float:
    filename = f"gotham.hydro_w.{row['sequence']}.bin"
    path = ARCHIVE_ROOT / row["input_relative"] / "node_00000000" / filename
    if not path.is_file():
        path = next((ARCHIVE_ROOT / row["input_relative"]).glob(f"node_*/{filename}"), None)
    if path is None or not path.is_file():
        raise RuntimeError(f"No representative cube found for {row}")
    with path.open("rb") as stream:
        for _ in range(8):
            line = stream.readline().decode("ascii")
            if line.strip().startswith("time="):
                return float(line.split("=", 1)[1])
    raise RuntimeError(f"Missing cube time in {path}")


def main() -> int:
    cubes, controls = read_catalog()
    mapped = 0
    unmapped = 0
    for manifest_path in sorted(MANIFEST_DIR.glob("res_*.tsv")):
        with manifest_path.open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if tuple(reader.fieldnames or ()) not in (OLD_FIELDS, NEW_FIELDS):
                raise RuntimeError(f"Unexpected header in {manifest_path}")
            rows = list(reader)

        for row in rows:
            key = (row["simulation"], row["phase"])
            cube_key = (key, row["sequence"])
            cube_time = cubes.get(cube_key)
            if cube_time is None:
                cube_time = read_cube_time(row)
            candidates = controls.get(key, [])
            if candidates:
                sequence, pdf_time = min(
                    candidates, key=lambda item: abs(item[1] - cube_time)
                )
                delta = abs(pdf_time - cube_time)
                if delta > time_tolerance(cube_time):
                    raise RuntimeError(
                        f"No matching PDF for {cube_key}: nearest delta={delta}"
                    )
                row["output_number"] = sequence
                row["output_time"] = repr(pdf_time)
                row["control_pdf_sequence"] = sequence
                mapped += 1
            else:
                row["output_number"] = row["sequence"]
                row["output_time"] = repr(cube_time)
                row["control_pdf_sequence"] = "none"
                unmapped += 1

        temporary = manifest_path.with_suffix(".tsv.partial")
        with temporary.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=NEW_FIELDS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(manifest_path)

    print(f"updated manifests: mapped={mapped}, without_archived_control={unmapped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
