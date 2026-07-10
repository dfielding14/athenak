#!/usr/bin/env python3
"""Generate traceable explanatory figures from existing GOTHAM PDF artifacts."""

from __future__ import annotations

import json
import math
import os
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


SCRIPT = Path(__file__).resolve()
DOCS_ROOT = SCRIPT.parents[1]
REPO_ROOT = SCRIPT.parents[3]
ANALYSIS_ROOT = Path(
    os.environ.get(
        "GOTHAM_ANALYSIS_ROOT",
        "/lustre/orion/ast207/proj-shared/gotham/analysis",
    )
)
ARCHIVE_PHASE = Path(
    "/lustre/orion/ast207/proj-shared/brent/gotham/res_8pc/phase2"
)
REBUILT = Path(
    "/lustre/orion/ast207/proj-shared/brent/gotham/"
    "pdf_rebuild_test/gpu_8pc_first8_original"
)
LOG_ROOT = Path(
    "/lustre/orion/ast207/proj-shared/brent/gotham/pdf_rebuild_test/logs"
)
GOLDEN_INPUT = Path(
    "/lustre/orion/ast207/proj-shared/brent/gotham/res_8pc/phase2/bin"
)

sys.path.insert(0, str(ANALYSIS_ROOT))
from gotham_analysis.readers.pdf import (  # noqa: E402
    _read_sparse_shard,
    read_pdf,
    read_pdf_header,
)


STREAMS = {
    "output29": "pdf_coord_abscostheta_coord_r",
    "output31": "pdf_coord_costheta_coord_r_temperature_vel_sph_r",
}


def output_directory(name):
    path = DOCS_ROOT / name / "figures"
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_figure(fig, directory, stem):
    for suffix in (".pdf", ".png"):
        fig.savefig(
            str(directory / (stem + suffix)),
            dpi=220 if suffix == ".png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def archived_first_n(stream, count=8):
    root = ARCHIVE_PHASE / "pdf" / stream
    header_path = root / "node_00000000" / "gotham.header.pdf"
    header = read_pdf_header(header_path)
    values = np.zeros(header["total_bins"], dtype=np.float64)
    times = set()
    source_files = []
    for index in range(count):
        path = root / ("node_{:08d}".format(index)) / "gotham.00062.pdf"
        time, indices, weights = _read_sparse_shard(path)
        times.add(time)
        np.add.at(values, indices.astype(np.intp), weights)
        source_files.append(str(path))
    if len(times) != 1:
        raise RuntimeError("archived sparse shard times differ: {!r}".format(times))
    return {
        "pdf": values.reshape(header["shape"]),
        "header": header,
        "time": next(iter(times)),
        "source_files": source_files,
        "header_path": str(header_path),
    }


def comparison_metrics(candidate, reference):
    reference_abs_sum = float(np.sum(np.abs(reference), dtype=np.float64))
    candidate_total = float(np.sum(candidate, dtype=np.float64))
    reference_total = float(np.sum(reference, dtype=np.float64))
    return {
        "e_l1": float(np.sum(np.abs(candidate - reference), dtype=np.float64))
        / reference_abs_sum,
        "e_total": abs(candidate_total - reference_total) / reference_abs_sum,
        "candidate_total": candidate_total,
        "reference_total": reference_total,
        "reference_abs_sum": reference_abs_sum,
    }


def load_controls():
    archived = {}
    rebuilt = {}
    for product, stream in STREAMS.items():
        archived[product] = archived_first_n(stream)
        payload = REBUILT / product / "gotham.00062.pdf"
        loaded = read_pdf(payload)
        rebuilt[product] = {
            "pdf": np.asarray(loaded["pdf"], dtype=np.float64),
            "time": float(loaded["time"]),
            "payload_path": str(payload),
            "header_path": str(REBUILT / product / "gotham.header.pdf"),
        }
    return archived, rebuilt


def safe_log_image(values):
    absolute = np.abs(values)
    positive = absolute[absolute > 0.0]
    floor = float(np.min(positive)) if positive.size else 1.0
    return np.maximum(absolute, floor), floor, float(np.max(absolute))


def figure_corruption_signature(archived, rebuilt, metrics):
    out = output_directory("corruption_root_cause")
    archive29 = archived["output29"]["pdf"]
    rebuilt29 = rebuilt["output29"]["pdf"]
    archive_abs, floor_a, max_a = safe_log_image(archive29)
    rebuilt_abs, floor_b, max_b = safe_log_image(rebuilt29)
    vmax = max(max_a, max_b)
    vmin = min(floor_a, floor_b)

    fig = plt.figure(figsize=(12.4, 8.4))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.78], hspace=0.33, wspace=0.3)
    axes = [fig.add_subplot(grid[0, index]) for index in range(3)]
    norm = LogNorm(vmin=vmin, vmax=vmax)
    images = [
        axes[0].imshow(archive_abs, origin="lower", aspect="auto", norm=norm, cmap="magma"),
        axes[1].imshow(rebuilt_abs, origin="lower", aspect="auto", norm=norm, cmap="magma"),
        axes[2].imshow(
            np.maximum(np.abs(rebuilt29 - archive29), vmin),
            origin="lower",
            aspect="auto",
            norm=norm,
            cmap="viridis",
        ),
    ]
    axes[0].set_title("Archived output29\nabsolute bin weight")
    axes[1].set_title("Rebuilt output29\nabsolute bin weight")
    axes[2].set_title("Absolute difference")
    for axis in axes:
        axis.set_xlabel("radius bin")
        axis.set_ylabel(r"$|\cos\theta|$ bin")
    fig.colorbar(images[1], ax=axes, fraction=0.025, pad=0.02, label="absolute weight")

    radial = fig.add_subplot(grid[1, :2])
    radial.plot(
        np.sum(archive29, axis=0, dtype=np.float64),
        color="#d95f02",
        lw=2.6,
        label="archived radial marginal",
    )
    radial.plot(
        np.sum(rebuilt29, axis=0, dtype=np.float64),
        color="#1b9e77",
        lw=1.7,
        ls="--",
        label="rebuilt radial marginal",
    )
    radial.axhline(0.0, color="0.5", lw=0.7)
    radial.set_xlabel("radius bin")
    radial.set_ylabel("signed net energy-flux weight")
    radial.set_title(
        r"Radial marginal agrees: $E_{\rm L1}=%.2e$"
        % metrics["output29_radial_marginal"]["e_l1"]
    )
    radial.legend(frameon=False, fontsize=9)

    angular = fig.add_subplot(grid[1, 2])
    x = np.arange(archive29.shape[0])
    archived_angle = np.sum(np.abs(archive29), axis=1, dtype=np.float64)
    rebuilt_angle = np.sum(np.abs(rebuilt29), axis=1, dtype=np.float64)
    angular.step(x, archived_angle, where="mid", color="#d95f02", lw=2.0, label="archive")
    angular.step(x, rebuilt_angle, where="mid", color="#1b9e77", lw=1.7, label="rebuild")
    angular.set_yscale("log")
    angular.set_xlabel(r"$|\cos\theta|$ bin")
    angular.set_ylabel("absolute weight")
    angular.set_title("Angular support does not agree")
    angular.legend(frameon=False, fontsize=9)
    fig.suptitle(
        "The diagnostic fingerprint: total and radial transport survive; angle does not",
        fontsize=14,
        y=0.99,
    )
    save_figure(fig, out, "output29_corruption_signature")


def figure_control_errors(metrics):
    out = output_directory("corruption_root_cause")
    names = [
        "output29 full\n(expected broken)",
        "output29 radial\nmarginal",
        "output31 full\n4D",
        "output31 geometry\nmarginal",
        "output31 temperature\nmarginal",
        "output31 radial velocity\nmarginal",
    ]
    keys = [
        "output29_full_expected_broken",
        "output29_radial_marginal",
        "output31_full_4d",
        "output31_geometry_costheta_radius",
        "output31_temperature_1d",
        "output31_radial_velocity_1d",
    ]
    values = [metrics[key]["e_l1"] for key in keys]
    colors = ["#d73027", "#1a9850", "#66bd63", "#66bd63", "#66bd63", "#66bd63"]
    fig, axis = plt.subplots(figsize=(10.8, 5.2))
    axis.bar(np.arange(len(values)), values, color=colors)
    axis.set_yscale("log")
    axis.set_xticks(np.arange(len(values)))
    axis.set_xticklabels(names, rotation=18, ha="right")
    axis.set_ylabel(r"normalized $L_1$ error")
    axis.set_title("Real eight-shard archived-control comparisons")
    axis.axhline(2.0e-3, color="0.35", ls="--", lw=1.0, label="full-validator output31 limit")
    axis.legend(frameon=False)
    axis.grid(axis="y", which="both", alpha=0.22)
    save_figure(fig, out, "archived_control_errors")


def parse_runtime_log(path):
    text = path.read_text(encoding="utf-8")
    match = re.search(
        r"processed shards=(\d+) blocks=(\d+) cells=(\d+) "
        r"payload_GiB=([0-9.eE+-]+) elapsed_s=([0-9.eE+-]+)",
        text,
    )
    if not match:
        raise RuntimeError("could not parse runtime log: {}".format(path))
    shards, blocks, cells = [int(value) for value in match.group(1, 2, 3)]
    gib, seconds = [float(value) for value in match.group(4, 5)]
    return {
        "path": str(path),
        "shards": shards,
        "blocks": blocks,
        "cells": cells,
        "payload_gib": gib,
        "elapsed_s": seconds,
        "gib_per_s": gib / seconds,
        "million_cells_per_s": cells / seconds / 1.0e6,
    }


def figure_runtime():
    out = output_directory("streaming_amr_reducer")
    runs = [
        parse_runtime_log(LOG_ROOT / "gpu_8pc_shard.3318736.out"),
        parse_runtime_log(LOG_ROOT / "gpu_8pc_8shard.3318742.out"),
    ]
    labels = ["1 K80 / 1 shard", "4 K80s / 8 shards"]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5))
    axes[0].bar(labels, [run["gib_per_s"] for run in runs], color=["#7570b3", "#1b9e77"])
    axes[0].set_ylabel("effective input throughput [GiB/s]")
    axes[0].set_title("Measured end-to-end throughput")
    axes[1].bar(
        labels,
        [run["million_cells_per_s"] for run in runs],
        color=["#7570b3", "#1b9e77"],
    )
    axes[1].set_ylabel("million cells/s")
    axes[1].set_title("Measured cell-processing rate")
    for axis in axes:
        axis.tick_params(axis="x", rotation=15)
        axis.grid(axis="y", alpha=0.22)
    fig.suptitle("Andes is a correctness platform, but the real-data path runs end to end")
    save_figure(fig, out, "andes_runtime_throughput")
    return runs


def figure_histogram_identities():
    out = output_directory("streaming_amr_reducer")
    manifest_path = REBUILT / "rebuild_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    groups = {}
    for product in manifest["products"]:
        groups.setdefault(product["weight"], []).append(product)
    rows = []
    for weight, products in sorted(groups.items()):
        if len(products) < 2:
            continue
        totals = np.asarray([product["sum"] for product in products], dtype=np.float64)
        scale = max(float(np.max(np.abs(totals))), np.finfo(np.float64).tiny)
        rows.append(
            {
                "weight": weight,
                "relative_spread": float((np.max(totals) - np.min(totals)) / scale),
                "products": [product["id"] for product in products],
            }
        )
    fig, axis = plt.subplots(figsize=(8.3, 4.6))
    axis.bar(
        [row["weight"] for row in rows],
        [row["relative_spread"] for row in rows],
        color="#1b9e77",
    )
    axis.set_yscale("log")
    axis.set_ylabel("relative spread in total histogram weight")
    axis.set_title("Products sharing a computed weight close to the same total")
    for index, row in enumerate(rows):
        if row["relative_spread"] == 0.0:
            axis.text(index, 9.5e-14, "exact 0", ha="center", va="bottom", fontsize=9)
    axis.grid(axis="y", which="both", alpha=0.22)
    save_figure(fig, out, "same_weight_total_closure")
    return {"manifest_path": str(manifest_path), "groups": rows}


def read_inventory():
    manifest_root = REPO_ROOT / "tools/gotham_pdf_rebuild/jobs/manifests"
    inventory = []
    for path in sorted(manifest_root.glob("*.tsv")):
        rows = path.read_text(encoding="utf-8").splitlines()[1:]
        mapped = sum(line.split("\t")[8] != "none" for line in rows)
        inventory.append(
            {
                "simulation": path.stem,
                "snapshots": len(rows),
                "mapped_controls": mapped,
                "without_control": len(rows) - mapped,
                "path": str(path),
            }
        )
    return inventory


def figure_inventory(inventory):
    out = output_directory("frontier_release_strategy")
    names = [row["simulation"].replace("_", "\n") for row in inventory]
    mapped = np.asarray([row["mapped_controls"] for row in inventory])
    missing = np.asarray([row["without_control"] for row in inventory])
    x = np.arange(len(names))
    fig, axis = plt.subplots(figsize=(10.2, 5.0))
    axis.bar(x, mapped, color="#1b9e77", label="has archived PDF control")
    axis.bar(x, missing, bottom=mapped, color="#e6ab02", label="no archived control")
    axis.set_xticks(x)
    axis.set_xticklabels(names)
    axis.set_ylabel("full-cube snapshots")
    axis.set_title("The production campaign is 161 separate reconstruction decisions")
    axis.legend(frameon=False)
    axis.grid(axis="y", alpha=0.22)
    save_figure(fig, out, "production_inventory")


def golden_input_inventory():
    pattern = "node_*/gotham.hydro_w.00028.bin"
    paths = sorted(GOLDEN_INPUT.glob(pattern))
    if not paths:
        raise RuntimeError("no golden input shards found")

    with paths[0].open("rb") as stream:
        location_size = variable_size = nvars = embedded_size = None
        while True:
            line = stream.readline()
            if not line:
                raise RuntimeError("truncated golden input preheader")
            if b"size of location=" in line:
                location_size = int(line.split(b"=", 1)[1])
            elif b"size of variable=" in line:
                variable_size = int(line.split(b"=", 1)[1])
            elif b"number of variables=" in line:
                nvars = int(line.split(b"=", 1)[1])
            elif b"header offset=" in line:
                embedded_size = int(line.split(b"=", 1)[1])
                stream.seek(embedded_size, 1)
                data_offset = stream.tell()
                break

    nx = 64
    record_size = 10 * 4 + 6 * location_size + nx ** 3 * nvars * variable_size
    blocks = 0
    total_bytes = 0
    for path in paths:
        size = path.stat().st_size
        payload = size - data_offset
        if payload < 0 or payload % record_size != 0:
            raise RuntimeError("inconsistent fixed-record shard: {}".format(path))
        blocks += payload // record_size
        total_bytes += size

    return {
        "source_glob": str(GOLDEN_INPUT / pattern),
        "representative_source": str(paths[0]),
        "shards": len(paths),
        "blocks": blocks,
        "cells": blocks * nx ** 3,
        "meshblock_cells": nx ** 3,
        "data_offset_bytes": data_offset,
        "record_size_bytes": record_size,
        "total_file_bytes": total_bytes,
    }


def figure_validation_ladder(runtime_runs, golden_inventory):
    out = output_directory("frontier_release_strategy")
    stages = [
        {"name": "one block\nall products", "cells": 64 ** 3, "status": "passed"},
        {"name": "one real\n8 pc shard", "cells": runtime_runs[0]["cells"], "status": "passed"},
        {"name": "eight real\n8 pc shards", "cells": runtime_runs[1]["cells"], "status": "passed"},
        {
            "name": "complete 1024-shard\nFrontier golden",
            "cells": golden_inventory["cells"],
            "status": "pending",
        },
    ]
    colors = ["#1b9e77" if stage["status"] == "passed" else "#e6ab02" for stage in stages]
    fig, axis = plt.subplots(figsize=(10.8, 4.8))
    axis.bar(
        np.arange(len(stages)),
        [stage["cells"] for stage in stages],
        color=colors,
    )
    axis.set_yscale("log")
    axis.set_xticks(np.arange(len(stages)))
    axis.set_xticklabels([stage["name"] for stage in stages])
    axis.set_ylabel("cells exercised")
    axis.set_title("Validation has crossed real-data and MPI/GPU scale; complete-snapshot checks remain")
    axis.grid(axis="y", which="both", alpha=0.22)
    axis.text(
        2,
        stages[2]["cells"] * 1.15,
        "passed",
        ha="center",
        color="#116530",
        fontweight="bold",
    )
    axis.text(
        3,
        stages[3]["cells"] * 1.15,
        "required release gate",
        ha="center",
        color="#9a6700",
        fontweight="bold",
    )
    save_figure(fig, out, "validation_ladder")
    return stages


def write_metrics(metrics):
    path = DOCS_ROOT / "figure_metrics.json"
    path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def main():
    archived, rebuilt = load_controls()
    metrics = {}

    archive29 = archived["output29"]["pdf"]
    rebuilt29 = rebuilt["output29"]["pdf"]
    archive31 = archived["output31"]["pdf"]
    rebuilt31 = rebuilt["output31"]["pdf"]
    metrics["comparisons"] = {
        "output29_full_expected_broken": comparison_metrics(rebuilt29, archive29),
        "output29_radial_marginal": comparison_metrics(
            np.sum(rebuilt29, axis=0, dtype=np.float64),
            np.sum(archive29, axis=0, dtype=np.float64),
        ),
        "output31_full_4d": comparison_metrics(rebuilt31, archive31),
        "output31_geometry_costheta_radius": comparison_metrics(
            np.sum(rebuilt31, axis=(2, 3), dtype=np.float64),
            np.sum(archive31, axis=(2, 3), dtype=np.float64),
        ),
        "output31_temperature_1d": comparison_metrics(
            np.sum(rebuilt31, axis=(0, 1, 3), dtype=np.float64),
            np.sum(archive31, axis=(0, 1, 3), dtype=np.float64),
        ),
        "output31_radial_velocity_1d": comparison_metrics(
            np.sum(rebuilt31, axis=(0, 1, 2), dtype=np.float64),
            np.sum(archive31, axis=(0, 1, 2), dtype=np.float64),
        ),
    }
    metrics["control_sources"] = {
        product: {
            "archived_header": archived[product]["header_path"],
            "archived_payloads": archived[product]["source_files"],
            "rebuilt_header": rebuilt[product]["header_path"],
            "rebuilt_payload": rebuilt[product]["payload_path"],
            "time": archived[product]["time"],
        }
        for product in STREAMS
    }
    figure_corruption_signature(archived, rebuilt, metrics["comparisons"])
    figure_control_errors(metrics["comparisons"])
    metrics["runtime_runs"] = figure_runtime()
    metrics["histogram_identity"] = figure_histogram_identities()
    metrics["inventory"] = read_inventory()
    figure_inventory(metrics["inventory"])
    metrics["golden_input_inventory"] = golden_input_inventory()
    metrics["validation_ladder"] = figure_validation_ladder(
        metrics["runtime_runs"], metrics["golden_input_inventory"]
    )
    metrics["generator"] = str(SCRIPT)
    metrics["analysis_reader"] = str(ANALYSIS_ROOT / "gotham_analysis/readers/pdf.py")
    metrics["notes"] = [
        "Only existing archived sparse PDF shards, rebuilt dense PDFs, logs, and manifests were read.",
        "The pending complete-golden cell count is an exact fixed-record inventory from live source-file sizes; it is not a completed reducer run.",
    ]
    path = write_metrics(metrics)
    print("wrote figures and {}".format(path))


if __name__ == "__main__":
    main()
