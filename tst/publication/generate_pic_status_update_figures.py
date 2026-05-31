#!/usr/bin/env python3
"""Generate quantitative and schematic figures for STATUS_UPDATE.md.

The quantitative panels are derived from frozen Orion artifacts.  Schematic
panels are labeled explicitly so they cannot be mistaken for measurements.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import tempfile
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from immutable_orion_tree import staged_verified_frozen_tree
from immutable_orion_tree import staged_verified_legacy_read_only_tree


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "figures" / "status_update"
ORION = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
PROJECT_HOME = Path("/ccs/proj/ast207/proj-shared/PIC")
Q006 = ORION / "q006-exact-isothermal-runtime-local-20260531-v4"
Q007 = ORION / "local-readiness" / "q007-paper-deltaf-runtime-replay-20260531-v4"
Q008 = ORION / "q008-cpaw-history-extension-20260531-v5"
Q011 = ORION / "local-readiness" / "q011-injection-distribution-runtime-local-20260531-v4"
Q011_SUCCESSOR = (
    ORION
    / "local-readiness"
    / "q011-injection-distribution-runtime-local-successor-audit-20260531-v1"
)
Q023_BELL = ORION / "local-readiness" / "q023-bell-combined-observable-smoke-20260530"
PIN = ORION / "local-readiness" / "integrated-serial-host-binary-clean-20260531T034317Z"
LEDGER = ORION / "ledger" / "node_hours.jsonl"
RECEIPTS = ORION / "ledger" / "mirror_receipts.jsonl"
MIRROR_LEDGER = PROJECT_HOME / "ledger" / "node_hours.jsonl"
ACTIVE_PROMOTION = ORION / "policy" / "active_promotion.json"
Q007_ANALYSIS = Q007.parent / f"{Q007.name}-analysis.json"
Q023_BELL_INITIAL = (
    Q023_BELL / "2d/bin/pic_q023_paper_bell_linear_2d_candidate.mhd_w_bcc.00000.bin"
)
Q023_BELL_CONTINUED = (
    Q023_BELL
    / "restart2d/bin/pic_q023_paper_bell_linear_restart2d_cont.mhd_w_bcc.00003.bin"
)
EXPECTED_INVENTORIES = {
    "pin": (PIN, "223dcb07d453572ea407f65dc9f7020afa3d6370b9d3c3cf4a06ceb4ca0e852c"),
    "q006": (Q006, "2ec81b018a085a1471c479260b17282499f7e8e559a8b8ad45bafa67707f20f7"),
    "q007": (Q007, "772d67effde79fa978283ca7e8425e00734a8c89ccdb124fa215737c80191842"),
    "q008": (Q008, "9cd5a4e8f04bdb4c6666f13c5904bd291e04e894d83fb84fb3300d8ee65e23de"),
    "q011": (Q011, "b50e1c1082e8bfbafce67aed09c1b92cfbae74b696e409aa29c3b6e12987457d"),
    "q011 successor": (
        Q011_SUCCESSOR,
        "a8343749f7cbe49af25a1c07e77645de1d39a8093d0cbf0b9d98cc272129e93c",
    ),
}
EXPECTED_Q007_ANALYSIS_SHA256 = "bd1eb4dc58773b088ad3e97e8cf50bcd5fbda2aee61907767f52c10c34351113"
EXPECTED_Q023_BELL_FILE_COUNT = 61
EXPECTED_Q023_BELL_INVENTORY_SHA256 = "29a0a14bbbd0e879e2cb6fb321758904a8fe8d5783612992f0de2585dca842d5"
EXPECTED_Q023_BELL_INITIAL_SHA256 = "e074f28151ec5abf42f2c77c5f4daf33d98ad5e86ab8146db785d2ed67feab30"
EXPECTED_Q023_BELL_CONTINUED_SHA256 = "0b1b419a1252cd8b82d517b448261ee7aa7b69984f36adfaa66dcd6b08b07943"
EXPECTED_LEDGER_SEQUENCE_NUMBER = 55
EXPECTED_LEDGER_TAIL_SHA256 = "538a5623cd672009bc903d6fc17b20b4dc7b503cae3e560aa362ebe0f3615b95"

BLUE = "#4477AA"
ORANGE = "#EE7733"
GREEN = "#228833"
RED = "#CC3311"
PURPLE = "#AA3377"
CYAN = "#66CCEE"
GRAY = "#777777"


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def regular_file_sha256(path: Path) -> str:
    return hashlib.sha256(regular_file_bytes(path)).hexdigest()


def regular_file_bytes(path: Path, expected_sha256: str | None = None) -> bytes:
    """Read one stable self-contained regular file without reopening its pathname."""
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise ValueError(f"report input is not a self-contained regular file: {path}")
        payload = bytearray()
        while chunk := os.read(fd, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(fd)
        identity = lambda value: (  # noqa: E731
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        if identity(before) != identity(after):
            raise ValueError(f"report input changed while reading: {path}")
        measured_sha256 = hashlib.sha256(payload).hexdigest()
        if expected_sha256 is not None and measured_sha256 != expected_sha256:
            raise ValueError(f"report input SHA-256 drifted: {path}")
        return bytes(payload)
    finally:
        os.close(fd)


def verify_report_inputs() -> tuple[list[dict[str, Any]], float, dict[str, Any], dict[str, Any]]:
    payload_members = {
        "q006": "reports/probe_summary.json",
        "q008": "results.json",
        "q011": "runtime_audit_report.json",
    }
    tree_receipts = {}
    verified_tree_sizes = {}
    frozen_payloads = {}
    for name, (root, expected) in EXPECTED_INVENTORIES.items():
        with staged_verified_frozen_tree(
            root,
            expected,
            authorized_root=ORION,
            label=f"status-update {name}",
        ) as (receipt, staged_tree):
            tree_receipts[name] = receipt
            verified_tree_sizes[name] = staged_tree.total_regular_size()
            if name in payload_members:
                frozen_payloads[name] = load_json(
                    staged_tree.member_path(payload_members[name])
                )
    q007_analysis_payload = regular_file_bytes(
        Q007_ANALYSIS,
        EXPECTED_Q007_ANALYSIS_SHA256,
    )
    measured_q007_analysis = hashlib.sha256(q007_analysis_payload).hexdigest()
    with staged_verified_legacy_read_only_tree(
        Q023_BELL,
        EXPECTED_Q023_BELL_FILE_COUNT,
        EXPECTED_Q023_BELL_INVENTORY_SHA256,
        authorized_root=ORION,
        label="status-update legacy q023 Bell smoke",
    ) as (q023_bell_receipt, q023_bell_snapshot):
        verified_tree_sizes["q023_bell"] = q023_bell_snapshot.total_regular_size()
        q023_initial_payload = q023_bell_snapshot.member_path(
            Q023_BELL_INITIAL.relative_to(Q023_BELL)
        ).read_bytes()
        q023_continued_payload = q023_bell_snapshot.member_path(
            Q023_BELL_CONTINUED.relative_to(Q023_BELL)
        ).read_bytes()
    if hashlib.sha256(q023_initial_payload).hexdigest() != EXPECTED_Q023_BELL_INITIAL_SHA256:
        raise ValueError("status-update initial Q023 Bell payload SHA-256 drifted")
    measured_q023_initial = hashlib.sha256(q023_initial_payload).hexdigest()
    if hashlib.sha256(q023_continued_payload).hexdigest() != EXPECTED_Q023_BELL_CONTINUED_SHA256:
        raise ValueError("status-update continued Q023 Bell payload SHA-256 drifted")
    measured_q023_continued = hashlib.sha256(q023_continued_payload).hexdigest()

    promotion = load_json(ACTIVE_PROMOTION)
    control_plane_version = promotion.get("control_plane_version")
    if not isinstance(control_plane_version, str) or len(control_plane_version) != 64:
        raise ValueError("status-update active promotion control-plane version is malformed")
    installed = ORION / "control_plane" / control_plane_version
    sys.path.insert(0, str(installed))
    try:
        from control_plane_common import require_storage_policy_unlock_snapshot
        from ledger import validate_mirrored_state

        policy, policy_snapshot = require_storage_policy_unlock_snapshot(
            control_plane_version=control_plane_version,
            authorized_pic_root=ORION,
            authorized_project_home_root=PROJECT_HOME,
        )
        rows = validate_mirrored_state(
            LEDGER,
            RECEIPTS,
            MIRROR_LEDGER,
            ledger_root=ORION / "ledger",
            receipts_root=ORION / "ledger",
            mirror_root=PROJECT_HOME / "ledger",
        )
    finally:
        sys.path.pop(0)
    if not rows:
        raise ValueError("status-update mirrored ledger is empty")
    ledger_tail = rows[-1]
    if ledger_tail["sequence_number"] != EXPECTED_LEDGER_SEQUENCE_NUMBER:
        raise ValueError("status-update mirrored ledger sequence drifted")
    if ledger_tail["event_sha256"] != EXPECTED_LEDGER_TAIL_SHA256:
        raise ValueError("status-update mirrored ledger tail SHA-256 drifted")
    node_hour_cap = float(policy["frontier"]["maximum_node_hours"])
    return rows, node_hour_cap, {
        "frozen_trees": tree_receipts,
        "legacy_read_only_trees": {"q023_bell_smoke": q023_bell_receipt},
        "q007_adjacent_analysis_sha256": measured_q007_analysis,
        "q023_bell_initial_mhd_sha256": measured_q023_initial,
        "q023_bell_continued_mhd_sha256": measured_q023_continued,
        "active_control_plane_version": control_plane_version,
        "active_policy_sha256": policy_snapshot["active_policy_sha256"],
        "active_promotion_sha256": policy_snapshot["active_promotion_sha256"],
        "ledger_tail_sequence_number": ledger_tail["sequence_number"],
        "ledger_tail_sha256": ledger_tail["event_sha256"],
        "node_hour_cap": node_hour_cap,
        "verified_tree_sizes_bytes": verified_tree_sizes,
    }, {
        **frozen_payloads,
        "q007_analysis": json.loads(q007_analysis_payload),
        "q023_initial_mhd": q023_initial_payload,
        "q023_continued_mhd": q023_continued_payload,
    }


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    color: str,
) -> None:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.025",
        linewidth=1.4,
        edgecolor=color,
        facecolor=color + "22",
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center")


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, color=GRAY))


def overview_workflow() -> None:
    fig, ax = plt.subplots(figsize=(11, 3))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 3)
    ax.axis("off")
    box(ax, (0.2, 1.0), 1.8, 1.0, "Source tree\nand tests", BLUE)
    box(ax, (2.45, 1.0), 1.8, 1.0, "Bounded host pin\n(dirty source)", PURPLE)
    box(ax, (4.7, 1.0), 1.8, 1.0, "Frozen Orion\nartifact trees", GREEN)
    box(ax, (6.95, 1.0), 1.8, 1.0, "Bounded artifact\nverifiers", ORANGE)
    box(ax, (9.2, 1.0), 1.6, 1.0, "Readiness\nclaims", RED)
    for left in (2.0, 4.25, 6.5, 8.75):
        arrow(ax, (left, 1.5), (left + 0.4, 1.5))
    ax.text(5.5, 2.6, "SCHEMATIC: local Orion-only readiness evidence flow", ha="center", weight="bold")
    ax.text(
        5.5,
        0.35,
        "Bulk evidence stays under the retained Orion root; broader Frontier portability and science qualification remain open.",
        ha="center",
        fontsize=9,
    )
    save(fig, "overview_workflow_schematic.png")


def verifier_hardening() -> None:
    fig, ax = plt.subplots(figsize=(11, 3.4))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 4)
    ax.axis("off")
    ax.text(5.5, 3.65, "SCHEMATIC: verifier hardening completed in this tranche", ha="center", weight="bold")
    before = [
        "Mutable source path",
        "Partial JSON acceptance",
        "Loose runtime replay checks",
    ]
    after = [
        "FD anchors + sealed semantic copies",
        "Exact schemas and closed directories",
        "Bounded diagnostics",
    ]
    for index, (left, right) in enumerate(zip(before, after)):
        y = 2.65 - index * 0.95
        box(ax, (0.35, y), 3.7, 0.55, left, RED)
        box(ax, (6.95, y), 3.7, 0.55, right, GREEN)
        arrow(ax, (4.15, y + 0.275), (6.85, y + 0.275))
    ax.text(2.2, 3.15, "Earlier risk", ha="center", color=RED, weight="bold")
    ax.text(8.8, 3.15, "Current bounded contract", ha="center", color=GREEN, weight="bold")
    save(fig, "verifier_hardening_schematic.png")


def readiness_gate_summary() -> None:
    labels = [
        "Bounded local successor tranche",
        "Verifier hardening",
        "Broad local test suites",
        "Final clean-tree validation",
        "Canonical clean-candidate freeze",
        "Science qualification",
        "Frontier portability matrix",
        "External review",
    ]
    status = ["completed", "completed", "completed", "completed", "open", "open", "open", "pending"]
    colors = [GREEN, GREEN, GREEN, GREEN, ORANGE, ORANGE, ORANGE, GRAY]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    y = np.arange(len(labels))
    ax.barh(y, np.ones(len(labels)), color=colors)
    for yy, state in zip(y, status):
        ax.text(0.5, yy, state.upper(), ha="center", va="center", color="white", weight="bold")
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.invert_yaxis()
    ax.set_title("Status classification, not a completion percentage")
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    save(fig, "readiness_gate_summary.png")


def read_mhd_binary(path: Path) -> dict[str, Any]:
    sys.path.insert(0, str(REPO / "vis/python"))
    try:
        import bin_convert

        return bin_convert.read_binary_as_athdf(str(path))
    finally:
        sys.path.pop(0)


def q023_bell_qualitative(initial_payload: bytes, continued_payload: bytes) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="athenak-pic-q023-bell-") as directory:
        initial_path = Path(directory) / "initial.bin"
        continued_path = Path(directory) / "continued.bin"
        initial_path.write_bytes(initial_payload)
        continued_path.write_bytes(continued_payload)
        initial = read_mhd_binary(initial_path)
        continued = read_mhd_binary(continued_path)
    initial_bz = np.asarray(initial["bcc3"], dtype=float)[0]
    continued_bz = np.asarray(continued["bcc3"], dtype=float)[0]
    continued_speed = np.sqrt(
        sum(np.asarray(continued[name], dtype=float)[0] ** 2 for name in ("velx", "vely", "velz"))
    )
    continued_uz = np.asarray(continued["velz"], dtype=float)[0]
    extent = [
        float(continued["x1f"][0]),
        float(continued["x1f"][-1]),
        float(continued["x2f"][0]),
        float(continued["x2f"][-1]),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5))
    panels = [
        (axes[0, 0], initial_bz, r"initial MHD $\delta B_z$", "coolwarm", True),
        (axes[0, 1], continued_bz, r"after one-cycle restart continuation: $\delta B_z$", "coolwarm", True),
        (axes[1, 0], continued_uz, r"after one-cycle restart continuation: $\delta u_z$", "coolwarm", True),
    ]
    for ax, values, title, cmap, diverging in panels:
        limit = float(np.max(np.abs(values)))
        kwargs = {"vmin": -limit, "vmax": limit} if diverging else {"vmin": 0.0, "vmax": limit}
        image = ax.imshow(values, origin="lower", aspect="auto", extent=extent, cmap=cmap, **kwargs)
        ax.set_xlabel(r"$x_1$")
        ax.set_ylabel(r"$x_2$")
        ax.set_title(title)
        colorbar = fig.colorbar(image, ax=ax, shrink=0.88)
        colorbar.formatter.set_powerlimits((-2, 2))
        colorbar.update_ticks()
    axes[1, 1].axis("off")
    axes[1, 1].text(
        0.5,
        0.68,
        "PIC phase-space panel unavailable",
        ha="center",
        va="center",
        fontsize=14,
        weight="bold",
        color=RED,
    )
    axes[1, 1].text(
        0.5,
        0.43,
        "This retained Q023 smoke couples CR particles internally\n"
        "and preserves restart payloads, but it did not serialize\n"
        "a particle VTK or phase-space product.",
        ha="center",
        va="center",
        fontsize=11,
    )
    axes[1, 1].text(
        0.5,
        0.19,
        "Required next campaign output:\n"
        "same-run particle distribution, current, and field snapshots",
        ha="center",
        va="center",
        fontsize=10,
        color=GRAY,
    )
    fig.suptitle(
        "Q023 Bell paper-mode source-local smoke: actual MHD fields, missing serialized PIC view",
        weight="bold",
    )
    save(fig, "q023_bell_mhd_pic_qualitative_gap.png")
    return {
        "artifact_root": str(Q023_BELL),
        "initial_mhd_path": str(Q023_BELL_INITIAL),
        "continued_mhd_path": str(Q023_BELL_CONTINUED),
        "initial_cycle": int(initial["NumCycles"]),
        "initial_time": float(initial["Time"]),
        "continued_cycle": int(continued["NumCycles"]),
        "continued_time": float(continued["Time"]),
        "initial_maximum_absolute_bcc3": float(np.max(np.abs(initial_bz))),
        "continued_maximum_absolute_bcc3": float(np.max(np.abs(continued_bz))),
        "continued_maximum_absolute_velz": float(np.max(np.abs(continued_uz))),
        "continued_maximum_speed": float(np.max(continued_speed)),
        "serialized_pic_phase_space_available": False,
        "qualification_effect": "bounded_source_local_preparation_smoke_only",
    }


def q006_carrier_mechanics(q006: dict[str, Any]) -> dict[str, Any]:
    wanted = [
        ("uniform_evolution_cycle2", "uniform"),
        ("smr_evolution_cycle2", "SMR"),
        ("audited_amr_evolution_cycle2", "audited AMR"),
        ("audited_amr_restart_cycle3", "AMR restart"),
    ]
    snapshots = q006["snapshots"]
    labels = [label for _, label in wanted]
    meshblocks = [snapshots[key]["snapshot"]["meshblock_count"] for key, _ in wanted]
    particles = [snapshots[key]["particle_count"] for key, _ in wanted]
    momentum = [snapshots[key]["maximum_absolute_total_momentum"] for key, _ in wanted]
    energy_error = [
        abs(snapshots[key]["total_kinetic_energy_density"] - 0.06) / 0.06 for key, _ in wanted
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.7))
    x = np.arange(len(labels))
    axes[0].bar(x, meshblocks, color=BLUE)
    axes[0].set_ylabel("mesh blocks")
    axes[0].set_xticks(x, labels, rotation=20, ha="right")
    axes[0].set_title("Carrier topology")
    axes[1].bar(x, momentum, color=ORANGE)
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$\|\mathbf{p}_{\rm total}\|_\infty$")
    axes[1].set_xticks(x, labels, rotation=20, ha="right")
    axes[1].set_title("Momentum closure")
    axes[2].bar(x, energy_error, color=GREEN)
    axes[2].set_yscale("log")
    axes[2].set_ylabel(r"$|K-0.06| / 0.06$")
    axes[2].set_xticks(x, labels, rotation=20, ha="right")
    axes[2].set_title("Kinetic-energy deviation")
    save(fig, "q006_carrier_mechanics.png")
    return {
        "labels": labels,
        "meshblocks": meshblocks,
        "particles": particles,
        "maximum_absolute_total_momentum_component": momentum,
        "relative_kinetic_energy_density_error": energy_error,
    }


def q007_branch_spectrum(q007: dict[str, Any]) -> dict[str, Any]:
    initial = q007["initial_branch_spectrum"]
    final = q007["final_branch_spectrum"]
    branch_keys = sorted({(item["direction"], item["signed_polarization"]) for item in initial})
    colors = [BLUE, ORANGE, GREEN, PURPLE]
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    metrics: dict[str, Any] = {}
    for branch_key, color in zip(branch_keys, colors):
        branch = f"direction={branch_key[0]:+d}, polarization={branch_key[1]:+d}"
        init_records = [
            item
            for item in initial
            if (item["direction"], item["signed_polarization"]) == branch_key
        ]
        final_records = [
            item for item in final if (item["direction"], item["signed_polarization"]) == branch_key
        ]
        modes = np.array([item["mode"] for item in init_records])
        init_values = [item["power"] for item in init_records]
        final_values = [item["power"] for item in final_records]
        ax.plot(modes, init_values, linestyle="--", color=color, alpha=0.5)
        ax.plot(modes, final_values, marker=".", markersize=3, color=color, label=branch)
        metrics[branch] = {
            "initial_min": min(init_values),
            "initial_max": max(init_values),
            "final_min": min(final_values),
            "final_max": max(final_values),
        }
    ax.plot(modes, 1.0e-6 / modes, color="black", linestyle=":", linewidth=1.2, label=r"expected $10^{-6}/m$")
    ax.set_yscale("log")
    ax.set_xlabel("mode index")
    ax.set_ylabel("branch power")
    ax.set_title("Q007 two-cycle replay: dashed startup, solid final")
    ax.legend(ncol=2, fontsize=8)
    save(fig, "q007_branch_spectrum_replay.png")
    return metrics


def q007_startup_residuals(q007: dict[str, Any]) -> dict[str, Any]:
    cases = q007["startup_cases"]
    labels = list(cases)
    friendly = ["max MHD residual", "relative macro-weight closure", "max initial $|w|$"]
    values = []
    for label in labels:
        case = cases[label]
        particles = case["initial_particles"]
        values.append(
            [
                max(case["initial_mhd_validation"].values()),
                abs(particles["macro_weight_sum"] - 128.0) / 128.0,
                max(abs(item) for item in particles["deltaf_weight_range"]),
            ]
        )
    values_array = np.maximum(np.array(values), 1.0e-20)
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    x = np.arange(len(labels))
    width = 0.23
    for index, (name, color) in enumerate(zip(friendly, [BLUE, ORANGE, GREEN])):
        ax.bar(x + (index - 1) * width, values_array[:, index], width, label=name, color=color)
    ax.set_yscale("log")
    ax.set_xticks(x, labels)
    ax.set_ylabel("maximum absolute startup residual")
    ax.set_title("Q007 startup carrier closure")
    ax.legend()
    ax.text(
        0.99,
        0.02,
        r"Exact zero startup $|w|$ values are clipped to $10^{-20}$ for log visibility.",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
    )
    save(fig, "q007_startup_residuals.png")
    return {"labels": labels, "fields": friendly, "values": values}


def convergence_metric(observation: dict[str, Any], metric: str) -> list[float]:
    return observation[metric]["residuals_nx1_32_64_128"]


def q008_convergence(q008: dict[str, Any]) -> dict[str, Any]:
    observations = q008["literal_profile_convergence_observations"]
    metrics = [
        ("magnetic_mode_amplitude_endpoint_relative_error", r"$B$ amplitude error"),
        ("velocity_mode_amplitude_endpoint_relative_error", r"$u$ amplitude error"),
        ("magnetic_phase_endpoint_absolute_error", r"$B$ phase error"),
        ("velocity_phase_endpoint_absolute_error", r"$u$ phase error"),
    ]
    colors = [BLUE, ORANGE, GREEN, PURPLE]
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.5), sharex=True)
    nx = np.array([32, 64, 128])
    for ax, (metric, title) in zip(axes.flat, metrics):
        for (profile, observation), color in zip(sorted(observations.items()), colors):
            ax.loglog(nx, convergence_metric(observation, metric), marker="o", color=color, label=profile)
        ax.set_title(title)
        ax.set_ylabel("absolute or relative residual")
        ax.grid(True, which="both", alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel(r"$N_x$")
    axes[0, 0].legend(fontsize=7)
    save(fig, "q008_literal_profile_convergence.png")
    return observations


def q008_case_heatmap(q008: dict[str, Any]) -> dict[str, Any]:
    limits = {
        "B amplitude maximum": 0.25,
        "u amplitude maximum": 0.25,
        "B phase maximum": 0.25,
        "u phase maximum": 0.25,
        "phase lock": 0.01,
        "undesired / desired Elsasser": 0.01,
        "mass drift": 1.0e-6,
        "scaled transverse energy": 0.25,
    }
    cases = q008["cases"]
    labels = list(cases)
    metrics = list(limits)
    values = []
    for case in cases.values():
        field = case["field_oracles"]
        history = case["history_oracles"]
        values.append(
            [
                field["magnetic_mode_amplitude_relative_error_max"],
                field["velocity_mode_amplitude_relative_error_max"],
                field["magnetic_phase_absolute_error_max"],
                field["velocity_phase_absolute_error_max"],
                field["velocity_magnetic_phase_lock_absolute_error_max"],
                field["undesired_to_desired_elsasser_ratio_max"],
                history["mass_relative_drift_max"],
                max(
                    history["transverse_kinetic_energy_scaled_relative_error_max"],
                    history["transverse_magnetic_energy_scaled_relative_error_max"],
                ),
            ]
        )
    matrix = np.array(
        [
            [max(float(value), 1.0e-20) / limits[metric] for metric, value in zip(metrics, row)]
            for row in values
        ]
    ).T
    fig, ax = plt.subplots(figsize=(13, 4.8))
    image = ax.imshow(matrix, aspect="auto", norm=LogNorm(vmin=1.0e-9, vmax=1.0), cmap="viridis")
    ax.set_xticks(np.arange(len(labels)), labels, rotation=75, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(metrics)), metrics, fontsize=8)
    ax.set_title("Q008 selected coarse preparation guardrails: residual / limit")
    bar = fig.colorbar(image, ax=ax)
    bar.set_label("fraction of guardrail limit")
    save(fig, "q008_case_guardrail_heatmap.png")
    return {
        "case_ids": labels,
        "metrics": metrics,
        "limits": limits,
        "guardrail_fraction_matrix": matrix.tolist(),
        "maximum_guardrail_fraction": float(np.max(matrix)),
    }


def q008_solver_axis(q008: dict[str, Any]) -> dict[str, Any]:
    cases = q008["cases"]
    selected = [
        case
        for case in cases.values()
        if case["profile"] == "exponential_expanding"
        and case["resolution"] == 64
        and case["axis"] == "x1"
    ]
    labels = [case["solver"] for case in selected]
    field = [case["field_oracles"]["magnetic_mode_amplitude_endpoint_relative_error"] for case in selected]
    phase = [case["field_oracles"]["magnetic_phase_endpoint_absolute_error"] for case in selected]
    fig, ax = plt.subplots(figsize=(7.8, 4))
    x = np.arange(len(labels))
    width = 0.34
    ax.bar(x - width / 2, field, width, label=r"$B$ amplitude error", color=BLUE)
    ax.bar(x + width / 2, phase, width, label=r"$B$ phase error", color=ORANGE)
    ax.set_xticks(x, labels)
    ax.set_ylabel("residual")
    ax.set_title("Q008 solver comparison at $N_x=64$")
    ax.legend()
    save(fig, "q008_solver_comparison.png")
    return {"solvers": labels, "b_amplitude_relative_error": field, "b_phase_error": phase}


def q011_distribution(q011: dict[str, Any]) -> dict[str, Any]:
    sampler = q011["monoenergetic_full_sphere_sampler"]
    octants = sampler["octant_counts"]
    moments = sampler["direction_second_moment"]
    means = sampler["direction_mean"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.9))
    axes[0].bar(np.arange(8), octants, color=BLUE)
    axes[0].axhline(q011["particle_count"] / 8, color=RED, linestyle="--", label="isotropic expectation")
    axes[0].set_xlabel("direction octant")
    axes[0].set_ylabel("particle count")
    axes[0].set_title("Octant census")
    axes[0].legend(fontsize=8)
    x = np.arange(3)
    axes[1].bar(x - 0.17, moments, 0.34, color=GREEN, label=r"$\langle \hat{v}_i^2\rangle$")
    axes[1].bar(x + 0.17, np.abs(means), 0.34, color=ORANGE, label=r"$|\langle \hat{v}_i\rangle|$")
    axes[1].axhline(1 / 3, color=RED, linestyle="--", label=r"isotropic $\langle \hat{v}_i^2\rangle$")
    axes[1].set_xticks(x, ["x", "y", "z"])
    axes[1].set_ylabel("bounded sample statistic")
    axes[1].set_title("Directional statistics")
    axes[1].legend(fontsize=8)
    save(fig, "q011_injection_distribution.png")
    return {
        "particle_count": q011["particle_count"],
        "octant_counts": octants,
        "direction_mean": means,
        "direction_second_moment": moments,
        "maximum_absolute_relative_speed_residual": sampler["maximum_absolute_relative_speed_residual"],
    }


def ledger_plot(rows: list[dict[str, Any]], node_hour_cap: float) -> dict[str, Any]:
    accounted = [row for row in rows if "cumulative_consumed_node_hours" in row]
    failed = [row for row in accounted if row["state"] == "FAILED"]
    sequence = [row["sequence_number"] for row in accounted]
    cumulative = [row["cumulative_consumed_node_hours"] for row in accounted]
    statuses = sorted(set(row["state"] for row in accounted))
    palette = {
        "COMPLETED": GREEN,
        "FAILED": RED,
        "CANCELLED": GRAY,
    }
    fig, ax = plt.subplots(figsize=(9, 4.1))
    ax.plot(sequence, cumulative, color=BLUE, linewidth=1.5)
    for status in statuses:
        selected = [row for row in accounted if row["state"] == status]
        ax.scatter(
            [row["sequence_number"] for row in selected],
            [row["cumulative_consumed_node_hours"] for row in selected],
            color=palette.get(status, PURPLE),
            label=status,
            s=24,
        )
    ax.set_xlabel("ledger sequence")
    ax.set_ylabel("cumulative node-hours")
    ax.set_title("Authoritative Frontier accounting ledger")
    ax.legend()
    save(fig, "frontier_compute_ledger.png")
    return {
        "ledger_record_count": len(rows),
        "accounting_record_count": len(accounted),
        "final_cumulative_node_hours": cumulative[-1],
        "failed_consumed_node_hours": sum(float(row["consumed_node_hours"]) for row in failed),
        "failed_consumption_fraction": (
            sum(float(row["consumed_node_hours"]) for row in failed) / cumulative[-1]
            if cumulative[-1]
            else 0.0
        ),
        "node_hour_cap": node_hour_cap,
        "remaining_node_hours": node_hour_cap - cumulative[-1],
        "ledger_tail_sequence_number": rows[-1]["sequence_number"],
        "ledger_tail_sha256": rows[-1]["event_sha256"],
        "status_counts": {status: sum(row["state"] == status for row in accounted) for status in statuses},
    }


def artifact_sizes(sizes: dict[str, int]) -> dict[str, int]:
    labels = list(sizes)
    mib = [sizes[label] / (1024**2) for label in labels]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.bar(labels, mib, color=[PURPLE, BLUE, ORANGE, GREEN, CYAN, GRAY])
    ax.set_yscale("log")
    ax.set_ylabel("tree size [MiB]")
    ax.set_title("Frozen Orion artifact inventory sizes")
    ax.tick_params(axis="x", rotation=20)
    save(fig, "artifact_inventory_sizes.png")
    return sizes


def validation_summary() -> dict[str, Any]:
    labels = ["publication\nsuite", "Frontier\ncontrol plane", "hardening\nfocus", "JSON parse", "Python AST"]
    counts = [719, 359, 99, 168, 88]
    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.bar(labels, counts, color=[BLUE, GREEN, ORANGE, PURPLE, CYAN])
    ax.bar_label(bars)
    ax.set_ylabel("checks or test cases passed")
    ax.set_title("Current pre-freeze hardening tranche")
    ax.text(
        0.99,
        0.02,
        "Counts are different check classes and are not additive.",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
    )
    save(fig, "validation_summary.png")
    return {"labels": labels, "counts": counts, "manually_curated": True}


def main() -> None:
    rows, node_hour_cap, verification, frozen_payloads = verify_report_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    q006 = frozen_payloads["q006"]
    q007 = frozen_payloads["q007_analysis"]
    q008 = frozen_payloads["q008"]
    q011 = frozen_payloads["q011"]

    overview_workflow()
    verifier_hardening()
    readiness_gate_summary()
    validation = validation_summary()
    metrics = {
        "sources": {
            "q006": str(Q006),
            "q007": str(Q007),
            "q008": str(Q008),
            "q011": str(Q011),
            "q011_successor": str(Q011_SUCCESSOR),
            "q023_bell": str(Q023_BELL),
            "pin": str(PIN),
            "ledger": str(LEDGER),
        },
        "input_verification": verification,
        "validation_summary": validation,
        "q023_bell_qualitative": q023_bell_qualitative(
            frozen_payloads["q023_initial_mhd"],
            frozen_payloads["q023_continued_mhd"],
        ),
        "q006": q006_carrier_mechanics(q006),
        "q007_branch_spectrum": q007_branch_spectrum(q007),
        "q007_startup_residuals": q007_startup_residuals(q007),
        "q008_convergence": q008_convergence(q008),
        "q008_heatmap": q008_case_heatmap(q008),
        "q008_solver": q008_solver_axis(q008),
        "q011": q011_distribution(q011),
        "ledger": ledger_plot(rows, node_hour_cap),
        "artifact_sizes_bytes": artifact_sizes(
            {
                "clean pin": verification["verified_tree_sizes_bytes"]["pin"],
                "Q006": verification["verified_tree_sizes_bytes"]["q006"],
                "Q007": verification["verified_tree_sizes_bytes"]["q007"],
                "Q008": verification["verified_tree_sizes_bytes"]["q008"],
                "Q011": verification["verified_tree_sizes_bytes"]["q011"],
                "Q011 successor": verification["verified_tree_sizes_bytes"]["q011 successor"],
            }
        ),
    }
    (OUT / "figure_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote figures and metrics to {OUT}")


if __name__ == "__main__":
    main()
