#!/usr/bin/env python3
"""Render a retained, nonqualifying review packet for Q011 pressure pilots."""

from __future__ import annotations

import argparse
import ctypes
import errno
import hashlib
import json
import os
from pathlib import Path
import stat
import uuid

import matplotlib

matplotlib.use("Agg")
from matplotlib import colors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

if __package__:
    from . import analyze_q011_section54_pressure_pilot as pilot
    from . import publish_q011_section54_pressure_pilot_bundle as publisher
else:
    import analyze_q011_section54_pressure_pilot as pilot
    import publish_q011_section54_pressure_pilot_bundle as publisher


PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
PUBLICATION_ROOT = PIC_ROOT / "publication"
DEFAULT_RECEIPT = PUBLICATION_ROOT / "q011_section54_pressure_pilot_bundle_receipt.json"
DEFAULT_OUTPUT = PUBLICATION_ROOT / "q011_section54_pressure_pilot_review_packet"
DEFAULT_OUTPUT_RECEIPT = (
    PUBLICATION_ROOT / "q011_section54_pressure_pilot_review_packet_receipt.json"
)
WATERMARK = "ENGINEERING CALIBRATION ONLY - NOT SUN & BAI REPRODUCTION"
_RENAME_NOREPLACE = 1


class PacketError(ValueError):
    """Raised when a review packet cannot be published safely."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PacketError(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )


def _direct_publication_path(value: str | Path, *, label: str) -> Path:
    path = Path(value)
    _require(path.is_absolute(), f"{label}: expected absolute path")
    lexical = Path(os.path.abspath(path))
    _require(lexical.parent == PUBLICATION_ROOT, f"{label}: expected direct publication child")
    return lexical


def _rename_no_replace(parent: Path, source_name: str, destination_name: str) -> None:
    descriptor = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        renameat2 = getattr(libc, "renameat2", None)
        _require(renameat2 is not None, "atomic no-replace rename is unavailable")
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        if (
            renameat2(
                descriptor,
                os.fsencode(source_name),
                descriptor,
                os.fsencode(destination_name),
                _RENAME_NOREPLACE,
            )
            != 0
        ):
            number = ctypes.get_errno()
            if number == errno.EEXIST:
                raise PacketError(f"publication target already exists: {destination_name}")
            raise OSError(number, os.strerror(number), destination_name)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _slice_xy(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim == 3:
        return array[array.shape[0] // 2]
    _require(array.ndim == 2, "expected a 2D or 3D retained field")
    return array


def _upstream_normalized(values: np.ndarray) -> np.ndarray:
    right = max(4, values.shape[1] // 10)
    reference = float(np.mean(values[:, -right:]))
    _require(np.isfinite(reference) and reference > 0.0, "upstream reference is invalid")
    return values / reference


def _terminal_products(
    bundle_root: Path, case: dict[str, object]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    snapshots = case["snapshots"]
    _require(type(snapshots) is list and bool(snapshots), "case snapshots are unavailable")
    snapshot = snapshots[-1]
    _require(type(snapshot) is dict, "terminal snapshot is malformed")
    time = float(snapshot["time"])
    ps_p0 = float(case["ps_p0"])
    _, mhd = pilot._binary_dataset(
        bundle_root, snapshot["mhd_w_bcc"], time, ps_p0, pilot._MHD_FIELDS
    )
    _, bmag = pilot._binary_dataset(
        bundle_root, snapshot["bmag"], time, ps_p0, pilot._SCALAR_BIN_FIELDS["bmag"]
    )
    _, jx = pilot._binary_dataset(
        bundle_root,
        snapshot["prtcl_jx"],
        time,
        ps_p0,
        pilot._SCALAR_BIN_FIELDS["prtcl_jx"],
    )
    dens_grid = mhd["dens"]
    return (
        np.asarray(dens_grid.x1_faces, dtype=float),
        np.asarray(dens_grid.x2_faces, dtype=float),
        _upstream_normalized(_slice_xy(dens_grid.values)),
        _upstream_normalized(_slice_xy(bmag["bmag"].values)),
        _slice_xy(jx["prtcl_jx"].values),
    )


def _qualitative_figure(bundle_root: Path, manifest: dict[str, object], output: Path) -> None:
    cases = manifest["cases"]
    _require(type(cases) is list and len(cases) == 4, "expected four pressure cases")
    products = [(case, _terminal_products(bundle_root, case)) for case in cases]
    max_current = max(
        float(np.percentile(np.abs(values[4]), 99.5)) for _, values in products
    )
    max_current = max(max_current, 1.0e-12)

    fig, axes = plt.subplots(4, 3, figsize=(14.2, 10.8), constrained_layout=True)
    for row, (case, (x1f, x2f, rho, bmag, jx)) in enumerate(products):
        title = f"{case['case_id']}: problem/ps_p0={case['ps_p0']}"
        bmin = max(0.5, float(np.min(bmag)))
        bmax = max(2.0, float(np.max(bmag)), bmin * (1.0 + 1.0e-12))
        density = axes[row, 0].pcolormesh(x1f, x2f, rho, shading="auto", cmap="viridis")
        magnetic = axes[row, 1].pcolormesh(
            x1f,
            x2f,
            np.maximum(bmag, 1.0e-12),
            shading="auto",
            cmap="magma",
            norm=colors.LogNorm(vmin=bmin, vmax=bmax),
        )
        current = axes[row, 2].pcolormesh(
            x1f,
            x2f,
            jx,
            shading="auto",
            cmap="coolwarm",
            norm=colors.TwoSlopeNorm(vmin=-max_current, vcenter=0.0, vmax=max_current),
        )
        axes[row, 0].set_ylabel(f"{title}\ny")
        for axis in axes[row]:
            axis.set_xlabel("x")
        fig.colorbar(density, ax=axes[row, 0], label="rho / upstream rho")
        fig.colorbar(magnetic, ax=axes[row, 1], label="|B| / upstream |B|")
        fig.colorbar(current, ax=axes[row, 2], label="PIC current jx")
    for axis, title in zip(
        axes[0], ("MHD density", "MHD magnetic magnitude", "PIC current diagnostic")
    ):
        axis.set_title(title)
    fig.suptitle(f"Q011 terminal pressure-pilot comparison\n{WATERMARK}", color="#8b0000")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def _profile_figure(analysis: dict[str, object], output: Path) -> None:
    records = [
        item
        for item in analysis["overlay_profiles"]
        if float(item["time_omega0_inverse"]) == 60.0
    ]
    _require(len(records) == 4, "expected four terminal overlay profiles")
    fig, axes = plt.subplots(2, 1, figsize=(10.4, 7.2), sharex=True, constrained_layout=True)
    for item in records:
        label = f"{item['case_id']} (ps_p0={item['ps_p0']})"
        x1 = np.asarray(item["x1_c_over_omega_pi"], dtype=float)
        axes[0].plot(x1, item["rho_y_average"], label=label)
        axes[1].plot(x1, item["bmag_y_average"], label=label)
    axes[0].set_ylabel("<rho>_y")
    axes[1].set_ylabel("<|B|>_y")
    axes[1].set_xlabel("x [c / omega_pi]")
    axes[0].legend(fontsize=8)
    axes[0].set_title(f"Q011 terminal pressure-pilot profiles\n{WATERMARK}", color="#8b0000")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def _summary_rows(analysis: dict[str, object]) -> list[dict[str, object]]:
    rows = []
    for case in analysis["case_summaries"]:
        telemetry = case["q017_telemetry"]
        rows.append(
            {
                "case_id": case["case_id"],
                "problem_ps_p0": case["ps_p0"],
                "terminal_particle_count": case["terminal_particles"]["particle_count"],
                "particle_efficiency": telemetry["load.particle_efficiency"],
                "zone_cycles_per_second": telemetry["throughput.zone_cycles_per_second"],
                "particle_updates_per_second": telemetry[
                    "throughput.particle_updates_per_second"
                ],
                "tracked_gpu_memory_high_water_bytes": telemetry[
                    "particle_memory.athenak_owned_tracked_kokkos_views."
                    "allocated_high_water_bytes_rank_max"
                ],
            }
        )
    return rows


def _markdown(rows: list[dict[str, object]]) -> bytes:
    records = [
        "# Q011 Section 5.4 pressure-pilot review packet",
        "",
        f"**{WATERMARK}**",
        "",
        "These four short runs are engineering calibration only. They do not select a pressure case automatically, authorize a qualifying launch, or establish a Sun and Bai reproduction claim.",
        "",
        "![Terminal MHD and PIC comparison](figures/terminal_mhd_pic_pressure_comparison.png)",
        "",
        "The qualitative panel compares retained terminal mesh density, magnetic-field magnitude, and the PIC current diagnostic for all four registered pressure choices.",
        "",
        "![Terminal profile overlays](figures/terminal_profile_overlays.png)",
        "",
        "The profile overlay compares terminal transverse averages without choosing a preferred case.",
        "",
        "| Case | `problem/ps_p0` | Terminal particles | Particle efficiency | Zone cycles/s | Particle updates/s | Tracked GPU memory high-water bytes |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        records.append(
            "| {case_id} | {problem_ps_p0} | {terminal_particle_count} | "
            "{particle_efficiency} | {zone_cycles_per_second} | "
            "{particle_updates_per_second} | {tracked_gpu_memory_high_water_bytes} |".format(
                **row
            )
        )
    records.extend(
        [
            "",
            "## Required human decision",
            "",
            "After reviewing the retained plots and metrics, record exactly one human-authored pressure-selection receipt for `ps_p0_1p00`, `ps_p0_0p05`, `ps_p0_0p10`, or `ps_p0_0p20`. The source-local validator intentionally does not rank or select cases.",
            "",
        ]
    )
    return ("\n".join(records)).encode("utf-8")


def render_packet(receipt_path: Path, output_path: Path, output_receipt: Path) -> dict[str, object]:
    receipt_path = _direct_publication_path(receipt_path, label="aggregate receipt")
    output_path = _direct_publication_path(output_path, label="packet output")
    output_receipt = _direct_publication_path(output_receipt, label="packet receipt")
    _require(not output_path.exists() and not output_receipt.exists(), "packet output already exists")
    publisher.verify_published_pressure_pilot_receipt(receipt_path)
    aggregate_receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    bundle_root = _direct_publication_path(
        aggregate_receipt["aggregate_bundle"]["path"], label="aggregate bundle"
    )
    analysis_path = _direct_publication_path(
        aggregate_receipt["aggregate_analysis"]["path"], label="aggregate analysis"
    )
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    manifest = pilot._manifest_schema(
        (bundle_root / pilot.MANIFEST_NAME).read_bytes()
    )

    staging = output_path.parent / f".{output_path.name}.staging-{uuid.uuid4()}"
    receipt_staging = output_receipt.parent / f".{output_receipt.name}.staging-{uuid.uuid4()}"
    staging.mkdir(mode=0o700)
    try:
        figures = staging / "figures"
        figures.mkdir(mode=0o700)
        qualitative = figures / "terminal_mhd_pic_pressure_comparison.png"
        profiles = figures / "terminal_profile_overlays.png"
        _qualitative_figure(bundle_root, manifest, qualitative)
        _profile_figure(analysis, profiles)
        rows = _summary_rows(analysis)
        metrics = staging / "pressure_review_metrics.json"
        metrics.write_bytes(
            _json_bytes(
                {
                    "schema_version": 1,
                    "record_type": "q011_section54_pressure_pilot_review_metrics",
                    "watermark": WATERMARK,
                    "qualification_effect": "none_no_selection_no_execution_authorization_no_science_claim",
                    "aggregate_receipt": {
                        "path": str(receipt_path),
                        "sha256": _sha256_path(receipt_path),
                    },
                    "cases": rows,
                }
            )
        )
        markdown = staging / "PRESSURE_REVIEW_PACKET.md"
        markdown.write_bytes(_markdown(rows))
        inventory_members = sorted(
            [
                qualitative.relative_to(staging),
                profiles.relative_to(staging),
                metrics.relative_to(staging),
                markdown.relative_to(staging),
            ],
            key=lambda path: path.as_posix(),
        )
        inventory = staging / "packet_inventory.json"
        inventory.write_bytes(
            _json_bytes(
                {
                    "schema_version": 1,
                    "record_type": "q011_section54_pressure_pilot_review_packet_inventory",
                    "members": [
                        {
                            "path": path.as_posix(),
                            "sha256": _sha256_path(staging / path),
                            "size": (staging / path).stat().st_size,
                        }
                        for path in inventory_members
                    ],
                }
            )
        )
        packet_receipt = {
            "schema_version": 1,
            "record_type": "q011_section54_pressure_pilot_review_packet_receipt",
            "watermark": WATERMARK,
            "qualification_effect": "none_no_selection_no_execution_authorization_no_science_claim",
            "aggregate_receipt": {
                "path": str(receipt_path),
                "sha256": _sha256_path(receipt_path),
            },
            "packet_root": str(output_path),
            "inventory_sha256": _sha256_path(inventory),
        }
        receipt_staging.write_bytes(_json_bytes(packet_receipt))
        for member in sorted(staging.rglob("*"), reverse=True):
            member.chmod(0o555 if member.is_dir() else 0o444)
        staging.chmod(0o555)
        receipt_staging.chmod(0o444)
        _rename_no_replace(output_path.parent, staging.name, output_path.name)
        _rename_no_replace(output_receipt.parent, receipt_staging.name, output_receipt.name)
        return packet_receipt
    except BaseException:
        if staging.exists():
            for member in sorted(staging.rglob("*"), reverse=True):
                member.chmod(0o700 if member.is_dir() else 0o600)
            staging.chmod(0o700)
            import shutil

            shutil.rmtree(staging)
        if receipt_staging.exists():
            receipt_staging.chmod(0o600)
            receipt_staging.unlink()
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aggregate-receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--output-receipt", type=Path, default=DEFAULT_OUTPUT_RECEIPT)
    args = parser.parse_args()
    result = render_packet(args.aggregate_receipt, args.output_path, args.output_receipt)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
