#!/usr/bin/env python3
"""Publish one immutable, nonqualifying Q011 pressure-pilot review packet."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile
from typing import Callable
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
QUALIFICATION_EFFECT = "none_no_selection_no_execution_authorization_no_science_claim"
INVENTORY_NAME = "packet_inventory.json"
PACKET_MEMBERS = (
    "PRESSURE_REVIEW_PACKET.md",
    "figures/terminal_mhd_pic_pressure_comparison.png",
    "figures/terminal_profile_overlays.png",
    "pressure_review_metrics.json",
)
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_FILE_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)


class PacketError(ValueError):
    """Raised when a review packet cannot be published safely."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PacketError(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )


def _direct_publication_path(
    value: str | Path, *, publication_root: Path, label: str
) -> Path:
    path = Path(value)
    _require(path.is_absolute(), f"{label}: expected absolute path")
    lexical = Path(os.path.abspath(path))
    _require(
        lexical.parent == publication_root and lexical.name not in {"", ".", ".."},
        f"{label}: expected direct publication child",
    )
    return lexical


def _read_readonly_at(parent_descriptor: int, name: str, label: str) -> bytes:
    with ImmutableReadonlyPublicationFile(parent_descriptor, name, label) as member:
        return member.read()


class ImmutableReadonlyPublicationFile:
    """Keep one direct publication member pinned while a consumer uses its bytes."""

    def __init__(self, parent_descriptor: int, name: str, label: str) -> None:
        _require("/" not in name, f"{label}: descriptor-relative read received nested path")
        self._parent_descriptor = parent_descriptor
        self._name = name
        self._label = label
        self._descriptor: int | None = None

    @property
    def descriptor(self) -> int:
        if self._descriptor is None:
            raise PacketError(f"{self._label} verifier is not open")
        return self._descriptor

    def __enter__(self) -> "ImmutableReadonlyPublicationFile":
        self._descriptor = os.open(
            self._name, _FILE_FLAGS, dir_fd=self._parent_descriptor
        )
        try:
            self.require_path_identity()
            _require(
                stat.S_ISREG(os.fstat(self.descriptor).st_mode)
                and not os.fstat(self.descriptor).st_mode & 0o222,
                f"{self._label} is not a read-only regular file",
            )
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *_: object) -> None:
        if self._descriptor is not None:
            os.close(self._descriptor)
            self._descriptor = None

    def require_path_identity(self) -> None:
        current = os.stat(
            self._name, dir_fd=self._parent_descriptor, follow_symlinks=False
        )
        opened = os.fstat(self.descriptor)
        _require(
            stat.S_ISREG(current.st_mode)
            and (current.st_dev, current.st_ino) == (opened.st_dev, opened.st_ino),
            f"{self._label} path changed while retained",
        )

    def read(self) -> bytes:
        self.require_path_identity()
        before = os.fstat(self.descriptor)
        os.lseek(self.descriptor, 0, os.SEEK_SET)
        payload = b""
        while True:
            chunk = os.read(self.descriptor, 1024 * 1024)
            if not chunk:
                break
            payload += chunk
        after = os.fstat(self.descriptor)
        _require(
            (
                before.st_dev,
                before.st_ino,
                before.st_size,
                before.st_mtime_ns,
            )
            == (
                after.st_dev,
                after.st_ino,
                after.st_size,
                after.st_mtime_ns,
            ),
            f"{self._label} changed while reading",
        )
        self.require_path_identity()
        return payload


class ImmutableReviewPacket:
    """Retain one no-follow packet descriptor while auditing its closure."""

    def __init__(
        self,
        root: Path,
        *,
        parent_fd: int,
        entry_name: str,
        inherited_root_fd: int | None = None,
    ) -> None:
        self.root = Path(os.path.abspath(root))
        self._parent_fd = parent_fd
        self._entry_name = entry_name
        self._inherited_root_fd = inherited_root_fd
        self._root_fd: int | None = None
        self._directory_identities: dict[str, tuple[int, int]] | None = None
        self._file_identities: dict[str, tuple[int, int]] | None = None

    @property
    def root_fd(self) -> int:
        if self._root_fd is None:
            raise PacketError("review-packet verifier is not open")
        return self._root_fd

    def __enter__(self) -> "ImmutableReviewPacket":
        self._root_fd = (
            os.open(self._entry_name, _DIRECTORY_FLAGS, dir_fd=self._parent_fd)
            if self._inherited_root_fd is None
            else os.dup(self._inherited_root_fd)
        )
        try:
            self.require_path_identity()
            _require(
                not os.fstat(self.root_fd).st_mode & 0o222,
                "review-packet root is writable",
            )
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *_: object) -> None:
        if self._root_fd is not None:
            os.close(self._root_fd)
            self._root_fd = None
        self._directory_identities = None
        self._file_identities = None

    def require_path_identity(self) -> None:
        publisher._require_same_directory_at(
            self._parent_fd,
            self._entry_name,
            self.root_fd,
            "review-packet root",
        )

    def _scan(
        self, directory_fd: int, prefix: tuple[str, ...] = ()
    ) -> tuple[list[str], dict[str, tuple[int, int]], dict[str, tuple[int, int]]]:
        paths = []
        directories = {}
        files = {}
        for name in sorted(os.listdir(directory_fd)):
            relative = "/".join((*prefix, name))
            metadata = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            if stat.S_ISREG(metadata.st_mode):
                _require(
                    not metadata.st_mode & 0o222,
                    f"review-packet member is writable: {relative}",
                )
                descriptor = os.open(name, _FILE_FLAGS, dir_fd=directory_fd)
                try:
                    opened = os.fstat(descriptor)
                    _require(
                        stat.S_ISREG(opened.st_mode)
                        and not opened.st_mode & 0o222
                        and (metadata.st_dev, metadata.st_ino)
                        == (opened.st_dev, opened.st_ino),
                        f"review-packet member changed during audit: {relative}",
                    )
                    files[relative] = (opened.st_dev, opened.st_ino)
                finally:
                    os.close(descriptor)
                paths.append(relative)
                continue
            _require(
                stat.S_ISDIR(metadata.st_mode) and not metadata.st_mode & 0o222,
                f"review-packet tree contains unsupported entry: {relative}",
            )
            descriptor = os.open(name, _DIRECTORY_FLAGS, dir_fd=directory_fd)
            try:
                opened = os.fstat(descriptor)
                _require(
                    (metadata.st_dev, metadata.st_ino)
                    == (opened.st_dev, opened.st_ino),
                    f"review-packet directory changed during audit: {relative}",
                )
                directories[relative] = (opened.st_dev, opened.st_ino)
                child_paths, child_directories, child_files = self._scan(
                    descriptor, (*prefix, name)
                )
                _require(
                    bool(child_paths),
                    f"review-packet tree contains empty directory: {relative}",
                )
                paths.extend(child_paths)
                directories.update(child_directories)
                files.update(child_files)
                current = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
                _require(
                    (current.st_dev, current.st_ino)
                    == (opened.st_dev, opened.st_ino),
                    f"review-packet directory changed during audit: {relative}",
                )
            finally:
                os.close(descriptor)
        return paths, directories, files

    def read(self, relative: str) -> bytes:
        self.require_path_identity()
        payload = publisher._ARTIFACT_HELPERS._read_at(
            self.root_fd,
            relative,
            self._directory_identities,
            self._file_identities,
        )
        self.require_path_identity()
        return payload

    def verify(self, expected_inventory_sha256: str) -> None:
        _require(
            publisher._SHA256_PATTERN.fullmatch(expected_inventory_sha256) is not None,
            "review-packet inventory SHA-256 is malformed",
        )
        paths, directories, files = self._scan(self.root_fd)
        _require(
            paths == sorted((*PACKET_MEMBERS, INVENTORY_NAME)),
            "review-packet tree closure drifted",
        )
        if self._directory_identities is None:
            self._directory_identities = directories
            self._file_identities = files
        else:
            _require(
                self._directory_identities == directories
                and self._file_identities == files,
                "review-packet member identity changed",
            )
        inventory_payload = self.read(INVENTORY_NAME)
        _require(
            _sha256(inventory_payload) == expected_inventory_sha256,
            "review-packet inventory SHA-256 drifted",
        )
        inventory = pilot._decode_json(inventory_payload, "review-packet inventory")
        _require(
            type(inventory) is dict
            and set(inventory) == {"schema_version", "record_type", "members"}
            and type(inventory["schema_version"]) is int
            and inventory["schema_version"] == 1
            and inventory["record_type"]
            == "q011_section54_pressure_pilot_review_packet_inventory"
            and type(inventory["members"]) is list,
            "review-packet inventory schema drifted",
        )
        members = inventory["members"]
        _require(
            [member.get("path") for member in members if type(member) is dict]
            == sorted(PACKET_MEMBERS),
            "review-packet inventory closure drifted",
        )
        for member in members:
            _require(
                type(member) is dict
                and set(member) == {"path", "sha256", "size"}
                and type(member["size"]) is int
                and member["size"] >= 0
                and type(member["sha256"]) is str
                and publisher._SHA256_PATTERN.fullmatch(member["sha256"]) is not None,
                "review-packet inventory member is malformed",
            )
            payload = self.read(member["path"])
            _require(
                len(payload) == member["size"] and _sha256(payload) == member["sha256"],
                f"review-packet member checksum drifted: {member['path']}",
            )
        self.require_path_identity()


def _verify_packet_at(
    root: Path,
    parent_descriptor: int,
    entry_name: str,
    root_descriptor: int,
    expected_inventory_sha256: str,
) -> None:
    with ImmutableReviewPacket(
        root,
        parent_fd=parent_descriptor,
        entry_name=entry_name,
        inherited_root_fd=root_descriptor,
    ) as packet:
        packet.verify(expected_inventory_sha256)
        packet.verify(expected_inventory_sha256)


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
    bundle_root: Path,
    case: dict[str, object],
    *,
    member_reader: Callable[[str], bytes],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    snapshots = case["snapshots"]
    _require(type(snapshots) is list and bool(snapshots), "case snapshots are unavailable")
    snapshot = snapshots[-1]
    _require(type(snapshot) is dict, "terminal snapshot is malformed")
    ps_p0 = float(case["ps_p0"])
    _, mhd = pilot._binary_dataset(
        bundle_root,
        snapshot["mhd_w_bcc"],
        ps_p0,
        pilot._MHD_FIELDS,
        member_reader=member_reader,
    )
    _, bmag = pilot._binary_dataset(
        bundle_root,
        snapshot["bmag"],
        ps_p0,
        pilot._SCALAR_BIN_FIELDS["bmag"],
        member_reader=member_reader,
    )
    _, jx = pilot._binary_dataset(
        bundle_root,
        snapshot["prtcl_jx"],
        ps_p0,
        pilot._SCALAR_BIN_FIELDS["prtcl_jx"],
        member_reader=member_reader,
    )
    dens_grid = mhd["dens"]
    return (
        np.asarray(dens_grid.x1_faces, dtype=float),
        np.asarray(dens_grid.x2_faces, dtype=float),
        _upstream_normalized(_slice_xy(dens_grid.values)),
        _upstream_normalized(_slice_xy(bmag["bmag"].values)),
        _slice_xy(jx["prtcl_jx"].values),
    )


def _qualitative_figure(
    bundle_root: Path,
    manifest: dict[str, object],
    output: Path,
    *,
    member_reader: Callable[[str], bytes],
) -> None:
    cases = manifest["cases"]
    _require(type(cases) is list and len(cases) == 4, "expected four pressure cases")
    products = [
        (case, _terminal_products(bundle_root, case, member_reader=member_reader))
        for case in cases
    ]
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


def _metrics(
    rows: list[dict[str, object]], aggregate_receipt: Path, receipt_payload: bytes
) -> bytes:
    return _json_bytes(
        {
            "schema_version": 1,
            "record_type": "q011_section54_pressure_pilot_review_metrics",
            "watermark": WATERMARK,
            "qualification_effect": QUALIFICATION_EFFECT,
            "aggregate_receipt": {
                "path": str(aggregate_receipt),
                "sha256": _sha256(receipt_payload),
            },
            "cases": rows,
        }
    )


def _inventory(payloads: dict[str, bytes]) -> bytes:
    return _json_bytes(
        {
            "schema_version": 1,
            "record_type": "q011_section54_pressure_pilot_review_packet_inventory",
            "members": [
                {"path": path, "sha256": _sha256(payloads[path]), "size": len(payloads[path])}
                for path in sorted(payloads)
            ],
        }
    )


def _verify_published_review_packet_receipt(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path = PIC_ROOT,
    allow_publication_guard: bool,
    require_publication_seal: bool = True,
) -> dict[str, str]:
    """Re-audit one durable packet receipt and its complete retained closure."""
    _pic_root, publication_root = publisher._publication_root(authorized_pic_root)
    acceptance_root = publisher._publication_acceptance_root(_pic_root)
    receipt = _direct_publication_path(
        receipt_path, publication_root=publication_root, label="packet receipt"
    )
    parent_descriptor = publisher._open_absolute_directory(publication_root)
    acceptance_descriptor = publisher._open_absolute_directory(acceptance_root)
    root_descriptor: int | None = None
    try:
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        publisher._require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "review-packet publication acceptance root",
        )
        if not allow_publication_guard:
            publisher._require_publication_guard_absent_at(
                parent_descriptor, receipt.name, "review-packet receipt"
            )
        receipt_identity = publisher._file_identity_at(
            parent_descriptor, receipt.name, "review-packet receipt"
        )
        receipt_payload = _read_readonly_at(
            parent_descriptor, receipt.name, "review-packet receipt"
        )
        if require_publication_seal:
            publisher._require_publication_seal_at(
                acceptance_descriptor,
                parent_descriptor,
                receipt.name,
                receipt_payload,
                receipt_identity,
                "review-packet receipt",
            )
        record = pilot._decode_json(receipt_payload, "review-packet receipt")
        _require(
            type(record) is dict
            and set(record)
            == {
                "schema_version",
                "record_type",
                "watermark",
                "qualification_effect",
                "consumption_rule",
                "publication_root_identity",
                "aggregate_receipt",
                "packet_root",
                "inventory_sha256",
                "source_bindings",
            }
            and type(record["schema_version"]) is int
            and record["schema_version"] == 1
            and record["record_type"]
            == "q011_section54_pressure_pilot_review_packet_receipt"
            and record["watermark"] == WATERMARK
            and record["qualification_effect"] == QUALIFICATION_EFFECT
            and record["consumption_rule"] == publisher.CONSUMPTION_RULE,
            "review-packet receipt schema drifted",
        )
        publisher._require_directory_identity(
            record["publication_root_identity"],
            parent_descriptor,
            "review-packet receipt publication root",
        )
        publisher._validate_retained_source_bindings(
            record["source_bindings"],
            require_verified_archive=publisher._is_production_pic_root(_pic_root),
        )
        aggregate = record["aggregate_receipt"]
        _require(
            type(aggregate) is dict
            and set(aggregate) == {"path", "sha256"}
            and type(aggregate["sha256"]) is str
            and publisher._SHA256_PATTERN.fullmatch(aggregate["sha256"]) is not None,
            "review-packet aggregate receipt binding is malformed",
        )
        aggregate_path = _direct_publication_path(
            aggregate["path"],
            publication_root=publication_root,
            label="review-packet aggregate receipt",
        )
        aggregate_payload = _read_readonly_at(
            parent_descriptor, aggregate_path.name, "retained aggregate receipt"
        )
        _require(
            _sha256(aggregate_payload) == aggregate["sha256"],
            "review-packet aggregate receipt drifted",
        )
        publisher.consume_published_pressure_pilot_bundle(
            aggregate_path, authorized_pic_root=authorized_pic_root
        )
        root = _direct_publication_path(
            record["packet_root"],
            publication_root=publication_root,
            label="review-packet root",
        )
        root_descriptor = os.open(root.name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor)
        _verify_packet_at(
            root,
            parent_descriptor,
            root.name,
            root_descriptor,
            record["inventory_sha256"],
        )
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        publisher._require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "review-packet publication acceptance root",
        )
        if not allow_publication_guard:
            publisher._require_publication_guard_absent_at(
                parent_descriptor, receipt.name, "review-packet receipt"
            )
        if require_publication_seal:
            publisher._require_publication_seal_at(
                acceptance_descriptor,
                parent_descriptor,
                receipt.name,
                receipt_payload,
                receipt_identity,
                "review-packet receipt",
            )
        return {
            "receipt_sha256": _sha256(receipt_payload),
            "inventory_sha256": record["inventory_sha256"],
            "aggregate_receipt_sha256": aggregate["sha256"],
        }
    finally:
        if root_descriptor is not None:
            os.close(root_descriptor)
        os.close(acceptance_descriptor)
        os.close(parent_descriptor)


def verify_published_review_packet_receipt(
    receipt_path: str | Path, *, authorized_pic_root: Path = PIC_ROOT
) -> dict[str, str]:
    """Re-audit one externally consumable packet receipt."""
    return _verify_published_review_packet_receipt(
        receipt_path,
        authorized_pic_root=authorized_pic_root,
        allow_publication_guard=False,
    )


def render_packet(
    receipt_path: Path,
    output_path: Path,
    output_receipt: Path,
    *,
    authorized_pic_root: Path = PIC_ROOT,
) -> dict[str, object]:
    """Render, freeze, verify, publish, and durably receipt one review packet."""
    _pic_root, publication_root = publisher._publication_root(authorized_pic_root)
    acceptance_root = publisher._publication_acceptance_root(_pic_root)
    receipt_path = _direct_publication_path(
        receipt_path, publication_root=publication_root, label="aggregate receipt"
    )
    output_path = _direct_publication_path(
        output_path, publication_root=publication_root, label="packet output"
    )
    output_receipt = _direct_publication_path(
        output_receipt, publication_root=publication_root, label="packet receipt"
    )
    parent_descriptor = publisher._open_absolute_directory(publication_root)
    acceptance_descriptor = publisher._open_absolute_directory(acceptance_root)
    staging = publication_root / f".{output_path.name}.staging-{uuid.uuid4()}"
    receipt_staging = publication_root / f".{output_receipt.name}.staging-{uuid.uuid4()}"
    staging_descriptor: int | None = None
    aggregate_receipt_file: ImmutableReadonlyPublicationFile | None = None
    aggregate_analysis_file: ImmutableReadonlyPublicationFile | None = None
    packet_renamed = False
    receipt_renamed = False
    guard_armed = False
    seal_committed = False
    receipt_identity: tuple[int, int] | None = None
    try:
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        publisher._require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "review-packet publication acceptance root",
        )
        publisher._require_absent_at(parent_descriptor, output_path.name, "packet output")
        publisher._require_absent_at(
            parent_descriptor, output_receipt.name, "packet receipt"
        )
        publisher._require_publication_guard_absent_at(
            parent_descriptor, output_receipt.name, "review-packet receipt"
        )
        publisher._require_publication_seal_absent_at(
            acceptance_descriptor, output_receipt.name, "review-packet receipt"
        )
        publisher.consume_published_pressure_pilot_bundle(
            receipt_path, authorized_pic_root=authorized_pic_root
        )
        aggregate_receipt_file = ImmutableReadonlyPublicationFile(
            parent_descriptor, receipt_path.name, "aggregate receipt"
        )
        aggregate_receipt_file.__enter__()
        aggregate_payload = aggregate_receipt_file.read()
        aggregate = pilot._decode_json(aggregate_payload, "aggregate receipt")
        bundle_root = _direct_publication_path(
            aggregate["aggregate_bundle"]["path"],
            publication_root=publication_root,
            label="aggregate bundle",
        )
        analysis_path = _direct_publication_path(
            aggregate["aggregate_analysis"]["path"],
            publication_root=publication_root,
            label="aggregate analysis",
        )
        aggregate_analysis_file = ImmutableReadonlyPublicationFile(
            parent_descriptor, analysis_path.name, "aggregate analysis"
        )
        aggregate_analysis_file.__enter__()
        analysis_payload = aggregate_analysis_file.read()
        _require(
            _sha256(analysis_payload) == aggregate["aggregate_analysis"]["sha256"],
            "aggregate analysis SHA-256 differs from aggregate receipt binding",
        )
        analysis = pilot._decode_json(analysis_payload, "aggregate analysis")
        bundle_descriptor = os.open(
            bundle_root.name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor
        )
        try:
            with publisher.ImmutablePressurePilotBundle(
                bundle_root,
                inherited_root_fd=bundle_descriptor,
                parent_fd=parent_descriptor,
                entry_name=bundle_root.name,
            ) as bundle:
                bundle.verify(aggregate["aggregate_bundle"]["manifest_sha256"])
                manifest = pilot._manifest_schema(bundle.read(pilot.MANIFEST_NAME))
                os.mkdir(staging.name, mode=0o700, dir_fd=parent_descriptor)
                staging_descriptor = os.open(
                    staging.name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor
                )
                with tempfile.TemporaryDirectory() as directory:
                    temporary = Path(directory)
                    qualitative = temporary / "terminal_mhd_pic_pressure_comparison.png"
                    profiles = temporary / "terminal_profile_overlays.png"
                    _qualitative_figure(
                        bundle_root,
                        manifest,
                        qualitative,
                        member_reader=bundle.read,
                    )
                    _profile_figure(analysis, profiles)
                    rows = _summary_rows(analysis)
                    payloads = {
                        "PRESSURE_REVIEW_PACKET.md": _markdown(rows),
                        "figures/terminal_mhd_pic_pressure_comparison.png": qualitative.read_bytes(),
                        "figures/terminal_profile_overlays.png": profiles.read_bytes(),
                        "pressure_review_metrics.json": _metrics(
                            rows, receipt_path, aggregate_payload
                        ),
                    }
                bundle.verify(aggregate["aggregate_bundle"]["manifest_sha256"])
        finally:
            os.close(bundle_descriptor)
        aggregate_analysis_file.require_path_identity()
        aggregate_receipt_file.require_path_identity()
        inventory_payload = _inventory(payloads)
        for relative, payload in {**payloads, INVENTORY_NAME: inventory_payload}.items():
            publisher._write_exclusive_below(staging_descriptor, relative, payload)
        publisher._freeze_anchored_tree(staging_descriptor)
        inventory_sha256 = _sha256(inventory_payload)
        _verify_packet_at(
            staging,
            parent_descriptor,
            staging.name,
            staging_descriptor,
            inventory_sha256,
        )
        packet_receipt = {
            "schema_version": 1,
            "record_type": "q011_section54_pressure_pilot_review_packet_receipt",
            "watermark": WATERMARK,
            "qualification_effect": QUALIFICATION_EFFECT,
            "consumption_rule": publisher.CONSUMPTION_RULE,
            "publication_root_identity": publisher._directory_identity(parent_descriptor),
            "aggregate_receipt": {
                "path": str(receipt_path),
                "sha256": _sha256(aggregate_payload),
            },
            "packet_root": str(output_path),
            "inventory_sha256": inventory_sha256,
            "source_bindings": publisher._source_bindings(
                require_verified_archive=publisher._is_production_pic_root(_pic_root)
            ),
        }
        packet_receipt_payload = _json_bytes(packet_receipt)
        publisher._write_exclusive_at(
            parent_descriptor, receipt_staging.name, packet_receipt_payload
        )
        receipt_identity = publisher._file_identity_at(
            parent_descriptor, receipt_staging.name, "review-packet staged receipt"
        )
        _require(
            aggregate_receipt_file.read() == aggregate_payload,
            "aggregate receipt changed while rendering review packet",
        )
        aggregate_analysis_file.require_path_identity()
        publisher.consume_published_pressure_pilot_bundle(
            receipt_path, authorized_pic_root=authorized_pic_root
        )
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        publisher._require_same_directory_at(
            parent_descriptor, staging.name, staging_descriptor, "review-packet staging tree"
        )
        publisher._require_absent_at(parent_descriptor, output_path.name, "packet output")
        publisher._rename_no_replace_at(
            parent_descriptor, staging.name, output_path.name
        )
        packet_renamed = True
        publisher._fsync_descriptor(parent_descriptor)
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        _verify_packet_at(
            output_path,
            parent_descriptor,
            output_path.name,
            staging_descriptor,
            inventory_sha256,
        )
        publisher._require_absent_at(
            parent_descriptor, output_receipt.name, "packet receipt"
        )
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        publisher._arm_publication_guard_at(parent_descriptor, output_receipt.name)
        guard_armed = True
        publisher._rename_no_replace_at(
            parent_descriptor, receipt_staging.name, output_receipt.name
        )
        receipt_renamed = True
        publisher._fsync_descriptor(parent_descriptor)
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        _verify_published_review_packet_receipt(
            output_receipt,
            authorized_pic_root=authorized_pic_root,
            allow_publication_guard=True,
            require_publication_seal=False,
        )
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        _require(receipt_identity is not None, "review-packet receipt identity is absent")
        publisher._require_same_file_at(
            parent_descriptor,
            output_receipt.name,
            receipt_identity,
            "canonical review-packet receipt",
        )
        if aggregate_analysis_file is not None:
            aggregate_analysis_file.__exit__(None, None, None)
            aggregate_analysis_file = None
        if aggregate_receipt_file is not None:
            aggregate_receipt_file.__exit__(None, None, None)
            aggregate_receipt_file = None
        publisher._disarm_publication_guard_at(parent_descriptor, output_receipt.name)
        guard_armed = False
        publisher._require_same_file_at(
            parent_descriptor,
            output_receipt.name,
            receipt_identity,
            "canonical review-packet receipt",
        )
        publisher._require_same_directory(
            publication_root, parent_descriptor, "review-packet publication root"
        )
        publisher._require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "review-packet publication acceptance root",
        )
        publisher._publish_publication_seal_at(
            acceptance_descriptor,
            parent_descriptor,
            output_receipt.name,
            packet_receipt_payload,
            receipt_identity,
        )
        seal_committed = True
        return packet_receipt
    except BaseException:
        if receipt_renamed and receipt_identity is not None:
            try:
                publisher._require_same_directory(
                    publication_root,
                    parent_descriptor,
                    "review-packet publication root",
                )
                publisher._require_same_directory(
                    acceptance_root,
                    acceptance_descriptor,
                    "review-packet publication acceptance root",
                )
                publisher._require_publication_seal_at(
                    acceptance_descriptor,
                    parent_descriptor,
                    output_receipt.name,
                    packet_receipt_payload,
                    receipt_identity,
                    "canonical review-packet receipt",
                )
                publisher._require_publication_guard_absent_at(
                    parent_descriptor,
                    output_receipt.name,
                    "canonical review-packet receipt",
                )
            except BaseException:
                pass
            else:
                seal_committed = True
                return packet_receipt
        rollback_error: BaseException | None = None
        if receipt_renamed:
            try:
                publisher._ensure_publication_guard_at(
                    parent_descriptor, output_receipt.name
                )
                guard_armed = True
            except BaseException as error:
                rollback_error = error
        for published, name, rollback, identity in (
            (
                receipt_renamed,
                output_receipt.name,
                publisher._rollback_published_file,
                receipt_identity,
            ),
            (
                packet_renamed,
                output_path.name,
                publisher._rollback_published_directory,
                staging_descriptor,
            ),
        ):
            if not published:
                continue
            try:
                _require(identity is not None, "review-packet rollback identity is absent")
                rollback(parent_descriptor, name, identity)
            except BaseException as error:
                rollback_error = rollback_error or error
        if not packet_renamed and staging_descriptor is not None:
            publisher._cleanup_anchored_tree(
                parent_descriptor,
                staging.name,
                staging_descriptor,
                "review-packet staging tree",
            )
        if not receipt_renamed:
            publisher._cleanup_anchored_file(
                parent_descriptor, receipt_staging.name, "review-packet staged receipt"
            )
        if guard_armed and rollback_error is None:
            try:
                publisher._disarm_publication_guard_at(
                    parent_descriptor, output_receipt.name
                )
                guard_armed = False
            except BaseException as error:
                rollback_error = error
        if rollback_error is not None:
            raise PacketError("cannot withdraw invalid review packet") from rollback_error
        raise
    finally:
        if aggregate_analysis_file is not None:
            aggregate_analysis_file.__exit__(None, None, None)
        if aggregate_receipt_file is not None:
            aggregate_receipt_file.__exit__(None, None, None)
        if staging_descriptor is not None:
            try:
                os.close(staging_descriptor)
            except OSError:
                if not seal_committed:
                    raise
        try:
            os.close(acceptance_descriptor)
            os.close(parent_descriptor)
        except OSError:
            if not seal_committed:
                raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aggregate-receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--output-receipt", type=Path, default=DEFAULT_OUTPUT_RECEIPT)
    parser.add_argument("--verify-published-receipt", type=Path)
    args = parser.parse_args()
    if args.verify_published_receipt is not None:
        verification = verify_published_review_packet_receipt(
            args.verify_published_receipt
        )
        print(json.dumps(verification, indent=2, sort_keys=True))
        return 0
    result = render_packet(args.aggregate_receipt, args.output_path, args.output_receipt)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
