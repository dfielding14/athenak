#!/usr/bin/env python3
"""Focused adversarial tests for the hardened Q011 applicability successor."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

try:
    from tst.publication import analyze_q011_section54_outputs as output_primitives
    from tst.publication import q011_section54_physical_applicability_successor_v1 as app
    from tst.publication import q011_section54_production_science_successor_v1 as science
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_physical_applicability_successor_v1 as app
    import q011_section54_production_science_successor_v1 as science


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_physical_applicability_successor_v1_2026-06-06.json"
)
DESIGN = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_physical_applicability_successor_v1_2026-06-06.md"
)
NX1 = 480
NX2 = 32
DX1 = 20.0
DX2 = 10.0
TIME = 500.0
SOURCE_COMMIT = app.TRUSTED_Q011_RUNTIME_SOURCE_COMMIT
ATTEMPT_ID = "q011-applicability-adversarial-fixture"
RAW_PATHS = {
    product: f"raw/{product}.00500.bin"
    for product in app.REQUIRED_RAW_PRODUCTS
}
RAW_PATHS["prtcl_all"] = "raw/prtcl_all.00500.part.vtk"


class _RetainedTemporaryDirectory(tempfile.TemporaryDirectory):
    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        self._finalizer.detach()


def _canonical(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _write_artifact(root: Path, role: str, relative: str, payload: bytes) -> dict[str, object]:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {
        "role": role,
        "path": relative,
        "sha256": _sha256_bytes(payload),
        "byte_count": len(payload),
    }


def _input_parameters(
    root_shape: tuple[int, int, int],
    block_shape: tuple[int, int, int],
    domain: tuple[float, float, float, float, float, float],
) -> dict[str, dict[str, str]]:
    return {
        "mesh": {
            "nx1": str(root_shape[0]),
            "nx2": str(root_shape[1]),
            "nx3": str(root_shape[2]),
            "nghost": "0",
            "x1min": str(domain[0]),
            "x1max": str(domain[1]),
            "x2min": str(domain[2]),
            "x2max": str(domain[3]),
            "x3min": str(domain[4]),
            "x3max": str(domain[5]),
        },
        "meshblock": {
            "nx1": str(block_shape[0]),
            "nx2": str(block_shape[1]),
            "nx3": str(block_shape[2]),
        },
    }


def _athena_binary_payload(dataset: output_primitives.AthenaBinaryDataset) -> bytes:
    parameter_header = "".join(
        f"<{block}>\n"
        + "".join(f"{key}={value}\n" for key, value in parameters.items())
        for block, parameters in dataset.input_parameters.items()
    ).encode("utf-8")
    payload = bytearray(
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        + f"  time={dataset.time}\n".encode("ascii")
        + f"  cycle={dataset.cycle}\n".encode("ascii")
        + f"  size of location={dataset.location_size}\n".encode("ascii")
        + f"  size of variable={dataset.variable_size}\n".encode("ascii")
        + f"  number of variables={len(dataset.variable_names)}\n".encode("ascii")
        + b"  variables:  "
        + b"  ".join(name.encode("ascii") for name in dataset.variable_names)
        + b"  \n"
        + f"  header offset={len(parameter_header)}\n".encode("ascii")
        + parameter_header
    )
    for block in dataset.blocks:
        payload.extend(struct.pack("<6i", *block.index_bounds))
        payload.extend(struct.pack("<4i", *block.logical_location, block.level))
        payload.extend(struct.pack("<6d", *block.geometry))
        values = np.concatenate(
            [
                np.asarray(block.fields[field], dtype="<f8").reshape(-1)
                for field in dataset.variable_names
            ]
        )
        payload.extend(values.tobytes())
    return bytes(payload)


def _block(
    shape: tuple[int, int, int],
    logical_location: tuple[int, int, int],
    level: int,
    geometry: tuple[float, float, float, float, float, float],
    fields: dict[str, np.ndarray],
) -> output_primitives.AthenaBinaryBlock:
    return output_primitives.AthenaBinaryBlock(
        index_bounds=(0, shape[0] - 1, 0, shape[1] - 1, 0, shape[2] - 1),
        logical_location=logical_location,
        level=level,
        geometry=geometry,
        fields=fields,
    )


def _spatial_fields(
    x: np.ndarray,
    y: np.ndarray,
    *,
    front: float,
    amplitude: float,
    high_k_fluctuation: bool,
    low_density_cell: bool = False,
) -> dict[str, np.ndarray]:
    front_index = int(np.searchsorted(x, front))
    density_x = np.full(x.size, 0.04 if high_k_fluctuation else 1.0)
    density_x[:front_index] = 4.0
    if front_index < x.size:
        density_x[front_index] = 3.8
    upstream = x > front
    density = np.broadcast_to(density_x[None, :], (y.size, x.size)).copy()
    if low_density_cell:
        density[0, -1] = 1.0e-4
    velx = np.broadcast_to(np.where(upstream, -30.0, 0.0)[None, :], density.shape).copy()
    zeros = np.zeros_like(density)
    wavelength = 40.0 if high_k_fluctuation else 400.0
    bcc2 = np.broadcast_to(
        amplitude * np.sin(2.0 * np.pi * x[None, :] / wavelength), density.shape
    ).copy()
    bcc3 = np.broadcast_to(
        amplitude * np.cos(2.0 * np.pi * y[:, None] / 160.0), density.shape
    ).copy()
    return {
        "dens": density[None, :, :],
        "velx": velx[None, :, :],
        "vely": zeros[None, :, :],
        "velz": zeros[None, :, :],
        "eint": np.full((1, y.size, x.size), 1.5),
        "bcc1": np.ones((1, y.size, x.size)),
        "bcc2": bcc2[None, :, :],
        "bcc3": bcc3[None, :, :],
    }


def _mhd_dataset(
    *,
    amplitude: float = 0.2,
    low_density_cell: bool = False,
    high_k_fluctuation: bool = False,
) -> output_primitives.AthenaBinaryDataset:
    domain = (0.0, NX1 * DX1, 0.0, NX2 * DX2, 0.0, 1.0)
    root_shape = (NX1, NX2, 1)
    x = (np.arange(NX1) + 0.5) * DX1
    y = (np.arange(NX2) + 0.5) * DX2
    fields = _spatial_fields(
        x,
        y,
        front=5000.0,
        amplitude=amplitude,
        high_k_fluctuation=high_k_fluctuation,
        low_density_cell=low_density_cell,
    )
    block = _block(
        (NX1, NX2, 1),
        (0, 0, 0),
        0,
        domain,
        fields,
    )
    return output_primitives.AthenaBinaryDataset(
        source=RAW_PATHS["mhd_w_bcc"],
        time=TIME,
        cycle=500,
        location_size=8,
        variable_size=8,
        variable_names=tuple(reversed(science.MHD_PRIMITIVE_FIELDS)),
        input_parameters=_input_parameters(root_shape, root_shape, domain),
        root_grid_shape=root_shape,
        meshblock_shape=root_shape,
        nghost=0,
        domain_bounds=domain,
        blocks=(block,),
    )


def _current_datasets(
    mhd: output_primitives.AthenaBinaryDataset,
    *,
    rho_q: float = 0.005,
    gas_frame_current: float = 0.02,
    rare_lambda_failure: bool = False,
) -> dict[str, output_primitives.AthenaBinaryDataset]:
    result: dict[str, output_primitives.AthenaBinaryDataset] = {}
    for product, field in science.CURRENT_PRODUCT_FIELDS.items():
        blocks = []
        for mhd_block in mhd.blocks:
            shape = mhd_block.fields["dens"].shape
            if product == "prtcl_rho":
                values = np.full(shape, rho_q)
            elif product == "prtcl_jx":
                values = rho_q * mhd_block.fields["velx"] + gas_frame_current
                if rare_lambda_failure:
                    values = values.copy()
                    values[0, -1, -1] = (
                        rho_q * mhd_block.fields["velx"][0, -1, -1] + 0.2
                    )
            else:
                values = np.zeros(shape)
            blocks.append(
                _block(
                    mhd.meshblock_shape,
                    mhd_block.logical_location,
                    mhd_block.level,
                    mhd_block.geometry,
                    {field: values},
                )
            )
        result[product] = output_primitives.AthenaBinaryDataset(
            **{
                **mhd.__dict__,
                "source": RAW_PATHS[product],
                "variable_names": (field,),
                "blocks": tuple(blocks),
            }
        )
    return result


def _mixed_level_mhd_dataset() -> output_primitives.AthenaBinaryDataset:
    root_shape = (480, 32, 1)
    block_shape = (240, 32, 1)
    domain = (0.0, 9600.0, 0.0, 320.0, 0.0, 1.0)
    leaves = (
        ((0, 0, 0), 1),
        ((0, 1, 0), 1),
        ((1, 0, 0), 1),
        ((1, 1, 0), 1),
        ((1, 0, 0), 0),
    )

    def geometry(location: tuple[int, int, int], level: int) -> tuple[float, ...]:
        result = []
        for axis, logical_index in enumerate(location):
            root_blocks = root_shape[axis] // block_shape[axis]
            logical_extent = root_blocks * (2**level if root_shape[axis] > 1 else 1)
            lower = domain[2 * axis]
            upper = domain[2 * axis + 1]
            width = (upper - lower) / logical_extent
            result.extend((lower + logical_index * width, lower + (logical_index + 1) * width))
        return tuple(result)

    blocks = []
    for location, level in leaves:
        bounds = geometry(location, level)
        x = np.linspace(bounds[0], bounds[1], block_shape[0], endpoint=False)
        x += 0.5 * (bounds[1] - bounds[0]) / block_shape[0]
        y = np.linspace(bounds[2], bounds[3], block_shape[1], endpoint=False)
        y += 0.5 * (bounds[3] - bounds[2]) / block_shape[1]
        blocks.append(
            _block(
                block_shape,
                location,
                level,
                bounds,
                _spatial_fields(
                    x,
                    y,
                    front=4500.0,
                    amplitude=0.2,
                    high_k_fluctuation=False,
                ),
            )
        )
    return output_primitives.AthenaBinaryDataset(
        source=RAW_PATHS["mhd_w_bcc"],
        time=TIME,
        cycle=500,
        location_size=8,
        variable_size=8,
        variable_names=tuple(reversed(science.MHD_PRIMITIVE_FIELDS)),
        input_parameters=_input_parameters(root_shape, block_shape, domain),
        root_grid_shape=root_shape,
        meshblock_shape=block_shape,
        nghost=0,
        domain_bounds=domain,
        blocks=tuple(blocks),
    )


def _particles(*, velocity: object = 10.0, x: float = 5500.0) -> dict[str, object]:
    count = app.PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES
    y = ((np.arange(count) % NX2) + 0.5) * DX2
    points = np.column_stack((np.full(count, x), y, np.full(count, 0.5)))
    values = np.asarray(velocity, dtype=np.float32)
    if values.ndim == 0:
        values = np.full(count, values, dtype=np.float32)
    velocities = np.column_stack((values, np.zeros(count), np.zeros(count))).astype(np.float64)
    return {
        "points": points,
        "cr_source": np.ones(count, dtype=np.int64),
        "birth_time": np.full(count, 45.0),
        "velocity": velocities,
        "macro_weight": np.ones(count),
    }


def _snapshot_particle_vtk(
    particles: dict[str, object], observed_time: float, cycle: int
) -> bytes:
    points = np.asarray(particles["points"], dtype=">f4")
    count = points.shape[0]
    integer_scalars = {
        "gid": np.zeros(count, dtype=">i4"),
        "ptag": np.arange(count, dtype=">i4"),
        "species": np.zeros(count, dtype=">i4"),
        "cr_source": np.asarray(particles["cr_source"], dtype=">i4"),
    }
    real_scalars = {
        "macro_weight": np.asarray(particles["macro_weight"], dtype=">f4"),
        "birth_time": np.asarray(particles["birth_time"], dtype=">f4"),
        "deltaf_f0": np.ones(count, dtype=">f4"),
        "deltaf_weight": np.zeros(count, dtype=">f4"),
    }
    velocities = np.asarray(particles["velocity"], dtype=">f4")
    payload = bytearray(
        (
            "# vtk DataFile Version 2.0\n"
            f"# AthenaK particle data at time= {observed_time}  nranks= 1  "
            f"cycle={cycle}  variables=prtcl_all\n"
            "BINARY\n"
            "DATASET UNSTRUCTURED_GRID\n"
            "\n"
            f"POINTS {count} float\n"
        ).encode("ascii")
    )
    payload.extend(points.tobytes())
    payload.extend(f"\n\nPOINT_DATA {count}\n".encode("ascii"))
    for name, values in integer_scalars.items():
        payload.extend(f"\nSCALARS {name} int\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(values.tobytes())
    for name, values in real_scalars.items():
        payload.extend(f"\nSCALARS {name} float\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(values.tobytes())
    payload.extend(b"\nVECTORS vel float\n")
    payload.extend(velocities.tobytes())
    return bytes(payload)


def _snapshot_evidence(
    root: Path,
    mhd: output_primitives.AthenaBinaryDataset,
    currents: dict[str, output_primitives.AthenaBinaryDataset],
    particles: dict[str, object],
    *,
    normalization: object | None = None,
    tamper_deck_after_binding: bool = False,
    ascii_raw_product: str | None = None,
) -> tuple[dict[str, object], dict[str, object]]:
    root.mkdir(parents=True, exist_ok=True)
    bindings = {
        "deck": _write_artifact(root, "deck", "bound/q011.in", b"<deck>\n"),
        "source_archive": _write_artifact(root, "source_archive", "bound/source.tar", b"source\n"),
        "executable": _write_artifact(root, "executable", "bound/athena", b"executable\n"),
    }
    source_manifest = {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.SOURCE_MANIFEST_RECORD_TYPE,
        "source_commit": SOURCE_COMMIT,
        "deck_sha256": bindings["deck"]["sha256"],
        "source_archive_sha256": bindings["source_archive"]["sha256"],
        "executable_sha256": bindings["executable"]["sha256"],
    }
    bindings["source_manifest"] = _write_artifact(
        root, "source_manifest", "bound/source-manifest.json", _canonical(source_manifest)
    )
    source_identity_receipt = {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.SOURCE_IDENTITY_RECEIPT_RECORD_TYPE,
        "authority": app.TRUSTED_IDENTITY_RECEIPT_AUTHORITY,
        "source_commit": SOURCE_COMMIT,
        "source_manifest_sha256": bindings["source_manifest"]["sha256"],
        "source_archive_sha256": bindings["source_archive"]["sha256"],
        "physical_escape_source_review_status": app.PHYSICAL_ESCAPE_SOURCE_REVIEW_STATUS,
    }
    bindings["source_identity_receipt"] = _write_artifact(
        root,
        "source_identity_receipt",
        "bound/source-identity-receipt.json",
        _canonical(source_identity_receipt),
    )
    executable_identity_receipt = {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.EXECUTABLE_IDENTITY_RECEIPT_RECORD_TYPE,
        "authority": app.TRUSTED_IDENTITY_RECEIPT_AUTHORITY,
        "source_commit": SOURCE_COMMIT,
        "source_manifest_sha256": bindings["source_manifest"]["sha256"],
        "source_identity_receipt_sha256": bindings["source_identity_receipt"]["sha256"],
        "executable_sha256": bindings["executable"]["sha256"],
        "build_identity_status": "digest_bound_candidate_only",
    }
    bindings["executable_identity_receipt"] = _write_artifact(
        root,
        "executable_identity_receipt",
        "bound/executable-identity-receipt.json",
        _canonical(executable_identity_receipt),
    )
    raw_payloads = {
        "mhd_w_bcc": _athena_binary_payload(mhd),
        **{
            product: _athena_binary_payload(currents[product])
            for product in science.CURRENT_PRODUCT_FIELDS
        },
        "prtcl_all": _snapshot_particle_vtk(particles, TIME, mhd.cycle),
    }
    if ascii_raw_product is not None:
        raw_payloads[ascii_raw_product] = b"arbitrary caller-controlled ASCII\n"
    raw = {
        product: _write_artifact(root, product, RAW_PATHS[product], raw_payloads[product])
        for product in app.REQUIRED_RAW_PRODUCTS
    }
    runtime_record = {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.NORMALIZATION_RECORD_TYPE,
        "source_commit": SOURCE_COMMIT,
        "normalization": (
            dict(app.EXACT_NORMALIZATION) if normalization is None else normalization
        ),
        "runtime_input_parameters_sha256": _sha256_bytes(_canonical(dict(mhd.input_parameters))),
        "deck_sha256": bindings["deck"]["sha256"],
        "source_manifest_sha256": bindings["source_manifest"]["sha256"],
        "source_identity_receipt_sha256": bindings["source_identity_receipt"]["sha256"],
        "source_archive_sha256": bindings["source_archive"]["sha256"],
        "executable_sha256": bindings["executable"]["sha256"],
        "executable_identity_receipt_sha256": bindings[
            "executable_identity_receipt"
        ]["sha256"],
    }
    runtime_binding = _write_artifact(
        root,
        "runtime_normalization_record",
        "bound/runtime-normalization.json",
        _canonical(runtime_record),
    )
    manifest = {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.SNAPSHOT_PROVENANCE_RECORD_TYPE,
        "attempt_id": ATTEMPT_ID,
        "source_commit": SOURCE_COMMIT,
        "executable_sha256": bindings["executable"]["sha256"],
        "runtime_normalization_sha256": runtime_binding["sha256"],
        "nominal_slot_time": TIME,
        "observed_committed_time": TIME,
        "cycle": mhd.cycle,
        "raw_products": raw,
        "decoded_product_sha256": {
            "mhd_w_bcc": app._dataset_sha256(mhd, "mhd_w_bcc"),
            **{
                product: app._dataset_sha256(currents[product], product)
                for product in science.CURRENT_PRODUCT_FIELDS
            },
            "prtcl_all": app._particle_payload_sha256(**particles),
        },
    }
    manifest_binding = _write_artifact(
        root, "snapshot_manifest", "bound/snapshot-manifest.json", _canonical(manifest)
    )
    if tamper_deck_after_binding:
        (root / bindings["deck"]["path"]).write_bytes(b"tampered deck\n")
    return (
        {
            "runtime_normalization_record": runtime_binding,
            **bindings,
        },
        {"snapshot_manifest": manifest_binding},
    )


def _snapshot(
    *,
    amplitude: float = 0.2,
    low_density_cell: bool = False,
    high_k_fluctuation: bool = False,
    rho_q: float = 0.005,
    gas_frame_current: float = 0.02,
    rare_lambda_failure: bool = False,
    particle_velocity: object = 10.0,
    particle_x: float = 5500.0,
    mixed_level: bool = False,
    normalization: object | None = None,
    tamper_deck_after_binding: bool = False,
    tamper_decoded_particle_after_manifest: bool = False,
    ascii_raw_product: str | None = None,
) -> app.ApplicabilitySnapshot:
    with _RetainedTemporaryDirectory() as temporary:
        root = Path(temporary)
        mhd = (
            _mixed_level_mhd_dataset()
            if mixed_level
            else _mhd_dataset(
                amplitude=amplitude,
                low_density_cell=low_density_cell,
                high_k_fluctuation=high_k_fluctuation,
            )
        )
        currents = _current_datasets(
            mhd,
            rho_q=rho_q,
            gas_frame_current=gas_frame_current,
            rare_lambda_failure=rare_lambda_failure,
        )
        particles = _particles(velocity=particle_velocity, x=particle_x)
        normalization_evidence, provenance = _snapshot_evidence(
            root,
            mhd,
            currents,
            particles,
            normalization=normalization,
            tamper_deck_after_binding=tamper_deck_after_binding,
            ascii_raw_product=ascii_raw_product,
        )
        if tamper_decoded_particle_after_manifest:
            particles["velocity"][0, 0] += 1.0
        return app.reduce_physical_applicability_snapshot(
            evidence_root=root,
            normalization_evidence=normalization_evidence,
            snapshot_provenance=provenance,
        )


def _state_vector(
    count: int = 0,
    weight: float = 0.0,
    energy: float = 0.0,
    momentum: list[float] | None = None,
) -> dict[str, object]:
    return {
        "particle_count": count,
        "macro_weight": weight,
        "kinetic_energy": energy,
        "momentum": [0.0, 0.0, 0.0] if momentum is None else momentum,
    }


def _face_state(
    reason_code: str,
    count: int = 0,
    weight: float = 0.0,
    energy: float = 0.0,
    momentum: list[float] | None = None,
) -> dict[str, object]:
    return {
        "reason_code": reason_code,
        **_state_vector(count, weight, energy, momentum),
    }


def _specific_energy_bounds(velocity: float = 10.0) -> tuple[float, float, float]:
    value = np.float32(velocity)
    restored = float(value)
    previous = float(np.nextafter(value, np.float32(-np.inf)))
    following = float(np.nextafter(value, np.float32(np.inf)))
    lower_velocity = restored - 0.5 * (restored - previous)
    upper_velocity = restored + 0.5 * (following - restored)

    def energy(speed: float) -> float:
        momentum_squared = speed * speed / (
            1.0 - speed * speed / science.PARTICLE_LIGHT_SPEED**2
        )
        return momentum_squared / (
            np.sqrt(1.0 + momentum_squared / science.PARTICLE_LIGHT_SPEED**2) + 1.0
        )

    return energy(lower_velocity), energy(restored), energy(upper_velocity)


def _checkpoint_velocities(count: int) -> np.ndarray:
    return np.resize(np.asarray([2.0, 4.0, 6.0, 8.0, 10.0], dtype=np.float32), count)


def _active_energy(count: int) -> float:
    return float(
        sum(_specific_energy_bounds(float(velocity))[1] for velocity in _checkpoint_velocities(count))
    )


def _active_energy_bins(count: int, edges: list[float]) -> tuple[list[int], list[float], list[float]]:
    specific = np.asarray(
        [
            _specific_energy_bounds(float(velocity))[1]
            for velocity in _checkpoint_velocities(count)
        ],
        dtype=np.float64,
    )
    indices = np.searchsorted(np.asarray(edges), specific, side="right") - 1
    indices = np.minimum(indices, len(edges) - 2)
    return (
        np.bincount(indices, minlength=len(edges) - 1).astype(int).tolist(),
        np.bincount(indices, weights=np.ones(count), minlength=len(edges) - 1).tolist(),
        np.bincount(indices, weights=specific, minlength=len(edges) - 1).tolist(),
    )


def _particle_checkpoint_vtk(
    count: int, observed_time: float, cycle: int, *, velocity: object | None = None
) -> bytes:
    points = np.zeros((count, 3), dtype=">f4")
    points[:, 0] = np.arange(count, dtype=np.float32)
    integer_scalars = {
        "gid": np.zeros(count, dtype=">i4"),
        "ptag": np.arange(count, dtype=">i4"),
        "species": np.zeros(count, dtype=">i4"),
        "cr_source": np.ones(count, dtype=">i4"),
    }
    real_scalars = {
        "macro_weight": np.ones(count, dtype=">f4"),
        "birth_time": np.full(count, app.STARTUP_REMOVAL_TIME, dtype=">f4"),
        "deltaf_f0": np.ones(count, dtype=">f4"),
        "deltaf_weight": np.zeros(count, dtype=">f4"),
    }
    velocities = np.zeros((count, 3), dtype=">f4")
    velocities[:, 0] = (
        _checkpoint_velocities(count)
        if velocity is None
        else np.asarray(velocity, dtype=np.float32)
    )
    payload = bytearray(
        (
            "# vtk DataFile Version 2.0\n"
            f"# AthenaK particle data at time= {observed_time}  nranks= 1  "
            f"cycle={cycle}  variables=prtcl_all\n"
            "BINARY\n"
            "DATASET UNSTRUCTURED_GRID\n"
            "\n"
            f"POINTS {count} float\n"
        ).encode("ascii")
    )
    payload.extend(points.tobytes())
    payload.extend(f"\n\nPOINT_DATA {count}\n".encode("ascii"))
    for name, values in integer_scalars.items():
        payload.extend(f"\nSCALARS {name} int\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(values.tobytes())
    for name, values in real_scalars.items():
        payload.extend(f"\nSCALARS {name} float\nLOOKUP_TABLE default\n".encode("ascii"))
        payload.extend(values.tobytes())
    payload.extend(b"\nVECTORS vel float\n")
    payload.extend(velocities.tobytes())
    return bytes(payload)


def _ps_escape_ledger(
    *,
    cycle: int,
    observed_time: float,
    active_count: int,
    active_mass: float,
    escaped_count: int = 0,
    escaped_mass: float = 0.0,
    escaped_energy: float = 0.0,
    escaped_momentum: list[float] | None = None,
    removed_count: int = 0,
    removed_mass: float = 0.0,
    initial_escape_count: int = 0,
) -> dict[str, object]:
    escaped_momentum = [0.0, 0.0, 0.0] if escaped_momentum is None else escaped_momentum
    injected_count = active_count + removed_count + escaped_count
    injected_mass = active_mass + removed_mass + escaped_mass
    return {
        "ps_cr_ledger_schema": app.PS_CR_LEDGER_SCHEMA,
        "ps_cr_ledger_complete": True,
        "ps_removed_excluded_early_cohort": True,
        "ps_escape_ledger_schema": app.PS_ESCAPE_LEDGER_SCHEMA,
        "ps_escape_ledger_complete": True,
        "ps_escape_audit_calls": app.PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE * cycle,
        "ps_escape_last_audit_time": observed_time,
        "ps_escaped_injected_cr_count_global": float(escaped_count),
        "ps_escaped_injected_cr_mass_global": escaped_mass,
        "ps_escaped_injected_cr_momentum_x1_global": escaped_momentum[0],
        "ps_escaped_injected_cr_momentum_x2_global": escaped_momentum[1],
        "ps_escaped_injected_cr_momentum_x3_global": escaped_momentum[2],
        "ps_escaped_injected_cr_energy_global": escaped_energy,
        "ps_escaped_initial_cr_count_global": float(initial_escape_count),
        "ps_injected_cr_count_global": float(injected_count),
        "ps_injected_cr_mass_global": injected_mass,
        "ps_injected_cr_momentum_x1_global": 0.0,
        "ps_injected_cr_momentum_x2_global": 0.0,
        "ps_injected_cr_momentum_x3_global": 0.0,
        "ps_injected_cr_energy_global": 1.0 if injected_count > 0 else 0.0,
        "ps_removed_cr_count_global": float(removed_count),
        "ps_removed_cr_mass_global": removed_mass,
        "ps_removed_cr_momentum_x1_global": 0.0,
        "ps_removed_cr_momentum_x2_global": 0.0,
        "ps_removed_cr_momentum_x3_global": 0.0,
        "ps_removed_cr_energy_global": 0.0,
    }


def _runtime_payload(snapshot: app.ApplicabilitySnapshot) -> dict[str, object]:
    record = snapshot.record
    identity_bindings = record["bound_normalization_evidence"]["bindings"]
    gates = record["gates"]
    particle = record["particle_gyroradius_containment"]
    exposure = record["particle_R_Lambda_exposure"]["populations"]["all_active"]
    metrics = {
        "R_maximum": gates["Q011-APP-R"]["observed_maximum"],
        "Lambda_maximum": gates["Q011-APP-LAMBDA"]["observed_maximum"],
        "S_delta_minimum_excluding_shock_transition": gates["Q011-APP-DI"][
            "observed_S_delta_minimum_excluding_shock_transition"
        ],
        "lambda_B_characteristic_over_local_di_maximum_minimum": gates["Q011-APP-DI"][
            "observed_lambda_B_characteristic_over_local_di_maximum"
        ],
        "sub_10di_magnetic_power_fraction_upper_bound_maximum": gates["Q011-APP-DI"][
            "observed_sub_10di_magnetic_power_fraction_upper_bound"
        ],
        "delta_B_rms_over_B0_minimum": gates["Q011-APP-DI"]["observed_delta_B_rms_over_B0"],
        "particle_rg_maximum_over_Ly": gates["Q011-APP-RG"]["observed_maximum_over_Ly"],
        "high_energy_tail_rg_maximum_over_Ly": gates["Q011-APP-RG"][
            "observed_high_energy_tail_maximum_over_Ly"
        ],
        "maximum_particle_specific_kinetic_energy": max(
            particle["maximum_specific_kinetic_energy"],
            _specific_energy_bounds()[2],
        ),
        "escaped_particle_rg_maximum_over_Ly": 0.0,
        "escaped_high_energy_tail_rg_maximum_over_Ly": 0.0,
        "escaped_particle_specific_kinetic_energy_maximum": 0.0,
    }
    inventory = []
    for index in range(int(app.EXPECTED_TERMINAL_TIME - app.STARTUP_REMOVAL_TIME) + 1):
        cycle = 45 + index
        if index == 0:
            start = 44.75
            end = app.STARTUP_REMOVAL_TIME
        else:
            start = app.STARTUP_REMOVAL_TIME + index - 1
            end = app.STARTUP_REMOVAL_TIME + index
        inventory.append(
            {
                "schema_version": app.SCHEMA_VERSION,
                "record_type": app.CYCLE_TELEMETRY_RECORD_TYPE,
                "attempt_id": record["snapshot_provenance"]["attempt_id"],
                "source_commit": record["snapshot_provenance"]["source_commit"],
                "executable_sha256": record["snapshot_provenance"]["executable_sha256"],
                "source_identity_receipt_sha256": identity_bindings[
                    "source_identity_receipt"
                ]["sha256"],
                "executable_identity_receipt_sha256": identity_bindings[
                    "executable_identity_receipt"
                ]["sha256"],
                "runtime_normalization_sha256": record["snapshot_provenance"][
                    "runtime_normalization_sha256"
                ],
                "ps_escape_accounting_source_commit": app.PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
                "cycle": cycle,
                "previous_committed_cycle": cycle - 1,
                "previous_committed_time": start,
                "start_time": start,
                "end_time": end,
                **metrics,
            }
        )
    count = particle["selected_particle_count"]
    active_energy = _active_energy(count)
    checkpoints = []
    for nominal in app.REQUIRED_PS_ESCAPE_CHECKPOINT_NOMINAL_TIMES:
        entry = next(item for item in inventory if item["end_time"] == nominal)
        cycle = entry["cycle"]
        active_count = count
        active_mass = float(count)
        checkpoints.append(
            {
                "nominal_checkpoint_time": nominal,
                "observed_committed_time": nominal,
                "cycle": cycle,
                "previous_committed_cycle": entry["previous_committed_cycle"],
                "previous_committed_time": entry["previous_committed_time"],
                "restart_artifact": None,
                "particle_checkpoint_artifact": None,
                "ps_escape_ledger": _ps_escape_ledger(
                    cycle=cycle,
                    observed_time=nominal,
                    active_count=active_count,
                    active_mass=active_mass,
                ),
                "active_injected_cr_count_global": float(active_count),
                "active_injected_cr_mass_global": active_mass,
                "active_injected_cr_kinetic_energy_global": active_energy,
                "escaped_injected_max_specific_kinetic_energy_global": 0.0,
                "escaped_injected_max_rg_over_Ly_global": 0.0,
            }
        )
    edges = [0.0, 5.0, 12.0, 25.0, 40.0, 60.0]
    counts, masses, energies = _active_energy_bins(count, edges)
    slope_cutoff_evidence = {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.SLOPE_CUTOFF_ESCAPE_RECORD_TYPE,
        "attempt_id": record["snapshot_provenance"]["attempt_id"],
        "source_commit": record["snapshot_provenance"]["source_commit"],
        "executable_sha256": record["snapshot_provenance"]["executable_sha256"],
        "source_identity_receipt_sha256": identity_bindings["source_identity_receipt"][
            "sha256"
        ],
        "executable_identity_receipt_sha256": identity_bindings[
            "executable_identity_receipt"
        ]["sha256"],
        "runtime_normalization_sha256": record["snapshot_provenance"][
            "runtime_normalization_sha256"
        ],
        "ps_escape_accounting_source_commit": app.PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
        "complete": True,
        "unavailable_reason": None,
        "high_energy_tail_threshold": 0.0,
        "energy_bin_edges": edges,
        "active_particle_count_by_bin": counts,
        "active_macro_weight_by_bin": masses,
        "active_kinetic_energy_by_bin": energies,
        "escaped_particle_count_by_bin": [0] * (len(edges) - 1),
        "escaped_macro_weight_by_bin": [0.0] * (len(edges) - 1),
        "escaped_kinetic_energy_by_bin": [0.0] * (len(edges) - 1),
        "active_particle_checkpoint_sha256": None,
        "escaped_particle_event_sha256": [],
    }
    return {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.RUNTIME_RECORD_TYPE,
        "successor_id": app.SUCCESSOR_ID,
        "qualification_effect": app.QUALIFICATION_EFFECT,
        "authorization": dict(app.AUTHORIZATION),
        "attempt_id": record["snapshot_provenance"]["attempt_id"],
        "source_commit": record["snapshot_provenance"]["source_commit"],
        "executable_sha256": record["snapshot_provenance"]["executable_sha256"],
        "source_identity_receipt_sha256": identity_bindings["source_identity_receipt"][
            "sha256"
        ],
        "executable_identity_receipt_sha256": identity_bindings[
            "executable_identity_receipt"
        ]["sha256"],
        "runtime_normalization_sha256": record["snapshot_provenance"][
            "runtime_normalization_sha256"
        ],
        "ps_escape_accounting_source_commit": app.PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
        "per_cycle_inventory": inventory,
        "ps_escape_checkpoints": checkpoints,
        "claim_escape_accounting": dict(app._CLAIM_ESCAPE_ACCOUNTING),
        "escaped_slope_cutoff_evidence": slope_cutoff_evidence,
        "escaped_particle_events": [],
        "particle_exposure": {
            "complete": True,
            "method": "every_particle_update_and_pre_destruction_boundary_event",
            "active_particle_updates_included": True,
            "pre_destruction_boundary_events_included": True,
            "escaped_particles_included": True,
            "escaped_particles_included_in_high_energy_and_gyroradius_extrema": True,
            "observation_count": count,
            "maximum_sampled_R": exposure["R"]["local_maximum"],
            "maximum_sampled_Lambda": exposure["Lambda"]["local_maximum"],
            "cumulative_macro_weighted_R_exceedance_fraction": exposure["R"][
                "macro_weighted"
            ]["exceedance_fraction"],
            "cumulative_CR_energy_weighted_R_exceedance_fraction": exposure["R"][
                "CR_kinetic_energy_weighted"
            ]["exceedance_fraction"],
            "cumulative_macro_weighted_Lambda_exceedance_fraction": exposure["Lambda"][
                "macro_weighted"
            ]["exceedance_fraction"],
            "cumulative_CR_energy_weighted_Lambda_exceedance_fraction": exposure[
                "Lambda"
            ]["CR_kinetic_energy_weighted"]["exceedance_fraction"],
        },
        "boundary_escape_ledger": {
            "complete": True,
            "scope": "all_Q011_shock_injected_particles_including_startup_removal",
            "nonperiodic_faces": {
                "inner_x1": _face_state(app.INNER_X1_ESCAPE_REASON),
                "outer_x1": _face_state(app.OUTER_X1_ESCAPE_REASON),
            },
            "periodic_faces": ["ix2", "ox2"],
            "accumulated_escaped": _state_vector(),
            "terminal_active": _state_vector(count, float(count), active_energy),
            "startup_removed": _state_vector(),
            "injected_particle_count": count,
            "injected_macro_weight": float(count),
            "escaped_particles_included_in_exposure": True,
            "escaped_particles_included_in_high_energy_and_gyroradius_extrema": True,
        },
    }


def _set_slope_cutoff_unavailable(
    payload: dict[str, object], reason: str = "escaped_binwise_tail_telemetry_unavailable"
) -> None:
    evidence = payload["escaped_slope_cutoff_evidence"]
    evidence["complete"] = False
    evidence["unavailable_reason"] = reason
    evidence["high_energy_tail_threshold"] = None
    for key in (
        "energy_bin_edges",
        "active_particle_count_by_bin",
        "active_macro_weight_by_bin",
        "active_kinetic_energy_by_bin",
        "escaped_particle_count_by_bin",
        "escaped_macro_weight_by_bin",
        "escaped_kinetic_energy_by_bin",
    ):
        evidence[key] = []


def _set_escape_events(
    payload: dict[str, object],
    *,
    count: int,
    macro_weight: float,
    kinetic_energy: float,
    rg_over_ly: float,
    momentum: list[float] | None = None,
    high_energy_tail_member: bool = False,
) -> None:
    momentum = [0.0, 0.0, 0.0] if momentum is None else momentum
    terminal = payload["per_cycle_inventory"][-1]
    event_weight = macro_weight / count
    specific_energy = kinetic_energy / macro_weight
    payload["escaped_particle_events"] = [
        {
            "schema_version": app.SCHEMA_VERSION,
            "record_type": app.ESCAPED_PARTICLE_EVENT_RECORD_TYPE,
            "attempt_id": payload["attempt_id"],
            "source_commit": payload["source_commit"],
            "executable_sha256": payload["executable_sha256"],
            "source_identity_receipt_sha256": payload["source_identity_receipt_sha256"],
            "executable_identity_receipt_sha256": payload[
                "executable_identity_receipt_sha256"
            ],
            "runtime_normalization_sha256": payload["runtime_normalization_sha256"],
            "ps_escape_accounting_source_commit": app.PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
            "event_index": index,
            "cycle": terminal["cycle"],
            "observed_committed_time": terminal["end_time"],
            "face": "outer_x1",
            "reason_code": app.OUTER_X1_ESCAPE_REASON,
            "particle_tag": 10000 + index,
            "macro_weight": event_weight,
            "specific_kinetic_energy": specific_energy,
            "macro_weighted_momentum": [component / count for component in momentum],
            "rg_over_Ly": rg_over_ly,
            "high_energy_tail_member": high_energy_tail_member,
        }
        for index in range(count)
    ]


def _configure_single_escape(
    payload: dict[str, object], *, specific_energy: float = 1.0, rg_over_ly: float = 0.01
) -> None:
    ledger = payload["boundary_escape_ledger"]
    active_energy = _active_energy(999)
    ledger["nonperiodic_faces"]["outer_x1"] = _face_state(
        app.OUTER_X1_ESCAPE_REASON, 1, 1.0, specific_energy
    )
    ledger["accumulated_escaped"] = _state_vector(1, 1.0, specific_energy)
    ledger["terminal_active"] = _state_vector(999, 999.0, active_energy)
    last_entry = payload["per_cycle_inventory"][-1]
    last_entry["maximum_particle_specific_kinetic_energy"] = max(
        last_entry["maximum_particle_specific_kinetic_energy"], specific_energy
    )
    last_entry["escaped_particle_specific_kinetic_energy_maximum"] = specific_energy
    last_entry["escaped_particle_rg_maximum_over_Ly"] = rg_over_ly
    last_entry["escaped_high_energy_tail_rg_maximum_over_Ly"] = rg_over_ly
    last_entry["particle_rg_maximum_over_Ly"] = max(
        last_entry["particle_rg_maximum_over_Ly"], rg_over_ly
    )
    last_entry["high_energy_tail_rg_maximum_over_Ly"] = max(
        last_entry["high_energy_tail_rg_maximum_over_Ly"], rg_over_ly
    )
    terminal = payload["ps_escape_checkpoints"][-1]
    terminal["active_injected_cr_count_global"] = 999.0
    terminal["active_injected_cr_mass_global"] = 999.0
    terminal["active_injected_cr_kinetic_energy_global"] = active_energy
    terminal["escaped_injected_max_specific_kinetic_energy_global"] = specific_energy
    terminal["escaped_injected_max_rg_over_Ly_global"] = rg_over_ly
    terminal["ps_escape_ledger"] = _ps_escape_ledger(
        cycle=terminal["cycle"],
        observed_time=terminal["observed_committed_time"],
        active_count=999,
        active_mass=999.0,
        escaped_count=1,
        escaped_mass=1.0,
        escaped_energy=specific_energy,
    )
    _set_escape_events(
        payload,
        count=1,
        macro_weight=1.0,
        kinetic_energy=specific_energy,
        rg_over_ly=rg_over_ly,
        high_energy_tail_member=True,
    )
    evidence = payload["escaped_slope_cutoff_evidence"]
    (
        evidence["active_particle_count_by_bin"],
        evidence["active_macro_weight_by_bin"],
        evidence["active_kinetic_energy_by_bin"],
    ) = _active_energy_bins(999, evidence["energy_bin_edges"])
    escaped_bin = int(
        np.searchsorted(evidence["energy_bin_edges"], specific_energy, side="right") - 1
    )
    escaped_bin = min(escaped_bin, len(evidence["energy_bin_edges"]) - 2)
    for key in (
        "escaped_particle_count_by_bin",
        "escaped_macro_weight_by_bin",
        "escaped_kinetic_energy_by_bin",
    ):
        evidence[key] = [0] * (len(evidence["energy_bin_edges"]) - 1)
    evidence["escaped_particle_count_by_bin"][escaped_bin] = 1
    evidence["escaped_macro_weight_by_bin"][escaped_bin] = 1.0
    evidence["escaped_kinetic_energy_by_bin"][escaped_bin] = specific_energy


def _shift_runtime_cycles(payload: dict[str, object], delta: int) -> None:
    for entry in payload["per_cycle_inventory"]:
        entry["cycle"] += delta
        entry["previous_committed_cycle"] += delta
    for checkpoint in payload["ps_escape_checkpoints"]:
        checkpoint["cycle"] += delta
        checkpoint["previous_committed_cycle"] += delta
        checkpoint["ps_escape_ledger"]["ps_escape_audit_calls"] = (
            app.PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE * checkpoint["cycle"]
        )


def _history(
    snapshots: list[app.ApplicabilitySnapshot],
    payload: dict[str, object],
    *,
    restart_ledger_overrides: dict[int, dict[str, object]] | None = None,
    particle_count_overrides: dict[int, int] | None = None,
    particle_payload_overrides: dict[int, bytes] | None = None,
    bind_cycle_telemetry: bool = True,
    bind_slope_cutoff_evidence: bool = True,
    bind_escape_events: bool = True,
) -> dict[str, object]:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        bound_payload = copy.deepcopy(payload)
        if bind_cycle_telemetry:
            bound_payload["per_cycle_inventory"] = [
                _write_artifact(
                    root,
                    "cycle_telemetry_record",
                    f"runtime/cycles/{index:05d}.json",
                    _canonical(entry),
                )
                for index, entry in enumerate(bound_payload["per_cycle_inventory"])
            ]
        for index, checkpoint in enumerate(bound_payload.get("ps_escape_checkpoints", [])):
            if checkpoint.get("restart_artifact") is None:
                restart_ledger = (
                    checkpoint["ps_escape_ledger"]
                    if restart_ledger_overrides is None
                    or index not in restart_ledger_overrides
                    else restart_ledger_overrides[index]
                )
                lines = ["<problem>"]
                for key, value in restart_ledger.items():
                    if type(value) is bool:
                        rendered = "true" if value else "false"
                    else:
                        rendered = str(value)
                    lines.append(f"{key} = {rendered}")
                lines.extend(("<par_end>", "restart-payload"))
                checkpoint["restart_artifact"] = _write_artifact(
                    root,
                    "ps_escape_restart_checkpoint",
                    f"checkpoints/q011.{index:02d}.rst",
                    ("\n".join(lines) + "\n").encode("latin1"),
                )
            if checkpoint.get("particle_checkpoint_artifact") is None:
                particle_count = int(checkpoint["active_injected_cr_count_global"])
                if particle_count_overrides is not None and index in particle_count_overrides:
                    particle_count = particle_count_overrides[index]
                checkpoint["particle_checkpoint_artifact"] = _write_artifact(
                    root,
                    "ps_escape_particle_checkpoint",
                    f"checkpoints/q011.{index:02d}.part.vtk",
                    (
                        particle_payload_overrides[index]
                        if particle_payload_overrides is not None
                        and index in particle_payload_overrides
                        else _particle_checkpoint_vtk(
                            particle_count,
                            float(checkpoint["observed_committed_time"]),
                            int(checkpoint["cycle"]),
                        )
                    ),
                )
        if bind_escape_events:
            bound_payload["escaped_particle_events"] = [
                _write_artifact(
                    root,
                    "escaped_particle_event",
                    f"runtime/escape-events/{index:05d}.json",
                    _canonical(event),
                )
                for index, event in enumerate(bound_payload["escaped_particle_events"])
            ]
        slope = bound_payload["escaped_slope_cutoff_evidence"]
        slope["active_particle_checkpoint_sha256"] = bound_payload[
            "ps_escape_checkpoints"
        ][-1]["particle_checkpoint_artifact"]["sha256"]
        if bind_escape_events:
            slope["escaped_particle_event_sha256"] = [
                binding["sha256"] for binding in bound_payload["escaped_particle_events"]
            ]
        if bind_slope_cutoff_evidence:
            bound_payload["escaped_slope_cutoff_evidence"] = _write_artifact(
                root,
                "slope_cutoff_escape_evidence",
                "runtime/slope-cutoff-escape.json",
                _canonical(slope),
            )
        binding = _write_artifact(
            root,
            "runtime_time_escape_record",
            "runtime/time-escape.json",
            _canonical(bound_payload),
        )
        return app.reduce_physical_applicability_history(
            snapshots,
            {"runtime_time_escape_record": binding},
            evidence_root=root,
        )


class Q011PhysicalApplicabilitySuccessorV1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.passing = _snapshot()

    def test_passing_snapshot_retains_full_bound_provenance_maps_regions_and_particles(self) -> None:
        record = self.passing.record
        self.assertTrue(record["snapshot_gate_pass_excluding_time_completeness"])
        self.assertEqual(set(self.passing.cell_maps), set(app.CELL_MAP_NAMES))
        self.assertEqual(set(record["regional_statistics"]), set(app.DETECTED_FRONT_REGIONS))
        self.assertEqual(
            set(record["snapshot_provenance"]["raw_products"]), set(app.REQUIRED_RAW_PRODUCTS)
        )
        for name, values in self.passing.cell_maps.items():
            self.assertFalse(values.flags.writeable)
            self.assertEqual(
                record["cell_map_contract"]["map_bindings"][name]["sha256"],
                app._array_sha256(values),
            )
        self.assertTrue(
            record["ion_scale_separation"]["precursor_magnetic_spectrum"]["actual_leaf_aware"]
        )
        self.assertFalse(record["authorization"]["claim_closure_authorized"])

    def test_normalization_is_bound_to_actual_artifacts_and_caller_dict_is_not_accepted(self) -> None:
        drifted = dict(app.EXACT_NORMALIZATION)
        drifted["selected_species_abs_q_over_m"] = 2.0
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "drifted"):
            _snapshot(normalization=drifted)
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "deck artifact .* drifted"):
            _snapshot(tamper_deck_after_binding=True)
        raw_only = _snapshot(tamper_decoded_particle_after_manifest=True)
        self.assertEqual(raw_only.record["gates"], self.passing.record["gates"])
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "not a trusted Athena binary"
        ):
            _snapshot(ascii_raw_product="mhd_w_bcc")
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "unexpected keyword"):
            app.reduce_physical_applicability_snapshot(  # type: ignore[call-arg]
                normalization=dict(app.EXACT_NORMALIZATION)
            )

    def test_history_rejects_stripped_snapshot_and_missing_full_evidence_sections(self) -> None:
        payload = _runtime_payload(self.passing)
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "complete ApplicabilitySnapshot"):
            _history([self.passing.record], payload)  # type: ignore[list-item]
        for section in (
            "snapshot_provenance",
            "regional_statistics",
            "ion_scale_separation",
            "particle_R_Lambda_exposure",
        ):
            record = copy.deepcopy(dict(self.passing.record))
            record.pop(section)
            stripped = app.ApplicabilitySnapshot(
                record=record,
                cell_maps=self.passing.cell_maps,
                evidence_root=self.passing.evidence_root,
            )
            with self.assertRaisesRegex(app.PhysicalApplicabilityError, "keys drifted"):
                _history([stripped], payload)

    def test_per_cell_R_and_Lambda_maxima_have_no_regional_dilution_waiver(self) -> None:
        r_failed = _snapshot(rho_q=0.02)
        self.assertFalse(r_failed.record["gates"]["Q011-APP-R"]["pass"])
        lambda_failed = _snapshot(rare_lambda_failure=True)
        gate = lambda_failed.record["gates"]["Q011-APP-LAMBDA"]
        regional = lambda_failed.record["regional_statistics"]["full_domain"]["Lambda"]
        self.assertFalse(gate["pass"])
        self.assertLess(regional["area_weighted"]["weighted_mean"], app.LAMBDA_MAXIMUM)
        for claim in (
            "Emax_claim",
            "high_energy_slope_or_cutoff_claim",
            "acceleration_rate_claim",
            "acceleration_efficiency_claim",
        ):
            self.assertIn(claim, lambda_failed.record["claim_rejections"])

    def test_DI_failures_reject_acceleration_claims_and_amplitude_floor_is_not_machine_tiny(self) -> None:
        for failed in (
            _snapshot(low_density_cell=True),
            _snapshot(high_k_fluctuation=True),
            _snapshot(amplitude=0.01),
        ):
            gate = failed.record["gates"]["Q011-APP-DI"]
            self.assertFalse(gate["pass"])
            for claim in (
                "Emax_claim",
                "high_energy_slope_or_cutoff_claim",
                "acceleration_rate_claim",
                "acceleration_efficiency_claim",
            ):
                self.assertIn(claim, failed.record["claim_rejections"])
        low_amplitude = _snapshot(amplitude=0.01)
        self.assertGreater(
            low_amplitude.record["gates"]["Q011-APP-DI"]["observed_delta_B_rms_over_B0"],
            np.finfo(np.float64).tiny,
        )
        self.assertLess(
            low_amplitude.record["gates"]["Q011-APP-DI"]["observed_delta_B_rms_over_B0"],
            app.DELTA_B_RMS_OVER_B0_MINIMUM,
        )

    def test_zero_fluctuation_power_fails_closed(self) -> None:
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "no resolvable fluctuation power"):
            _snapshot(amplitude=0.0)

    def test_high_energy_and_global_maxima_defeat_q999_bulk_dilution(self) -> None:
        velocities = np.full(app.PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES, 10.0)
        velocities[-1] = 100.0
        failed = _snapshot(particle_velocity=velocities)
        particle = failed.record["particle_gyroradius_containment"]
        self.assertLess(particle["macro_q999_over_Ly"], app.RG_MAXIMUM_OVER_LY_MAXIMUM)
        self.assertGreater(particle["maximum_over_Ly"], app.RG_MAXIMUM_OVER_LY_MAXIMUM)
        self.assertGreater(
            particle["energy_fraction_with_rg_above_Ly_over_4"],
            app.RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM,
        )
        self.assertFalse(failed.record["gates"]["Q011-APP-RG"]["pass"])
        self.assertIn("Emax_claim", failed.record["claim_rejections"])

    def test_mixed_level_AMR_uses_actual_leaf_aware_estimator_not_finest_composite_FFT(self) -> None:
        snapshot = _snapshot(mixed_level=True)
        spectrum = snapshot.record["ion_scale_separation"]["precursor_magnetic_spectrum"]
        self.assertEqual(spectrum["source_levels_present"], [0, 1])
        self.assertEqual(spectrum["analysis_source_level"], 0)
        self.assertEqual(spectrum["restriction_factor"], 2)
        self.assertTrue(spectrum["actual_leaf_aware"])
        self.assertTrue(spectrum["finest_composite_FFT_rejected"])

    def test_bound_complete_actual_per_cycle_inventory_passes_without_authorizing(self) -> None:
        history = _history([self.passing], _runtime_payload(self.passing))
        self.assertTrue(history["all_physical_applicability_gates_pass"])
        coverage = history["runtime_time_escape_evidence"]["cycle_coverage"]
        self.assertEqual(coverage["covered_cycle_count"], 1156)
        self.assertEqual(coverage["first_interval_crossing_start_time"], 44.75)
        self.assertEqual(coverage["first_interval_crossing_end_time"], 45.0)
        self.assertEqual(coverage["post_startup_removal_start_time"], 45.0)
        self.assertTrue(
            coverage["continuous_coverage_begins_exactly_at_startup_removal_time"]
        )
        self.assertTrue(coverage["continuous_coverage_includes_first_startup_crossing"])
        self.assertEqual(coverage["previous_committed_time_before_startup_crossing"], 44.75)
        self.assertEqual(coverage["terminal_time"], 1200.0)
        ps_escape = history["runtime_time_escape_evidence"]["ps_escape_accounting"]
        self.assertEqual(
            ps_escape["source_commit"], app.PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT
        )
        self.assertEqual(
            [item["nominal_checkpoint_time"] for item in ps_escape["checkpoints"]],
            list(app.REQUIRED_PS_ESCAPE_CHECKPOINT_NOMINAL_TIMES),
        )
        self.assertEqual(ps_escape["checkpoints"][0]["nominal_checkpoint_time"], 100.0)
        terminal = ps_escape["terminal"]
        self.assertEqual(
            terminal["ps_escape_ledger"]["ps_escape_audit_calls"],
            app.PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE * terminal["cycle"],
        )
        self.assertEqual(
            terminal["ps_escape_ledger"]["ps_escaped_initial_cr_count_global"], 0.0
        )
        self.assertFalse(history["authorization"]["claim_closure_authorized"])

    def test_ps_escape_checkpoint_schema_cadence_audit_time_and_restart_binding_fail_closed(self) -> None:
        missing = _runtime_payload(self.passing)
        missing["ps_escape_checkpoints"].pop()
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "cadence is incomplete"):
            _history([self.passing], missing)

        incomplete = _runtime_payload(self.passing)
        incomplete["ps_escape_checkpoints"][5]["ps_escape_ledger"][
            "ps_escape_ledger_complete"
        ] = False
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "is incomplete"):
            _history([self.passing], incomplete)

        schema = _runtime_payload(self.passing)
        schema["ps_escape_checkpoints"][5]["ps_escape_ledger"][
            "ps_escape_ledger_schema"
        ] = 2
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "schema-1 ps_escape"):
            _history([self.passing], schema)

        audit = _runtime_payload(self.passing)
        audit["ps_escape_checkpoints"][5]["ps_escape_ledger"]["ps_escape_audit_calls"] += 1
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "audit-call coverage"):
            _history([self.passing], audit)

        time_drift = _runtime_payload(self.passing)
        time_drift["ps_escape_checkpoints"][5]["ps_escape_ledger"][
            "ps_escape_last_audit_time"
        ] += 0.5
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "last audit time"):
            _history([self.passing], time_drift)

        forged = _runtime_payload(self.passing)
        forged_restart = copy.deepcopy(forged["ps_escape_checkpoints"][5]["ps_escape_ledger"])
        forged_restart["ps_escape_audit_calls"] += 2
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "disagrees with bound restart"):
            _history([self.passing], forged, restart_ledger_overrides={5: forged_restart})

        unbound = _runtime_payload(self.passing)
        unbound["ps_escape_checkpoints"][5]["restart_artifact"] = {
            "role": "ps_escape_restart_checkpoint",
            "path": "missing/unbound.rst",
            "sha256": "0" * 64,
            "byte_count": 1,
        }
        with self.assertRaises(app.PhysicalApplicabilityError):
            _history([self.passing], unbound)

        missing_particle = _runtime_payload(self.passing)
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "active count disagrees with bound particle"
        ):
            _history([self.passing], missing_particle, particle_count_overrides={5: 999})

        particle_time = _runtime_payload(self.passing)
        checkpoint = particle_time["ps_escape_checkpoints"][5]
        payload = _particle_checkpoint_vtk(
            int(checkpoint["active_injected_cr_count_global"]),
            float(checkpoint["observed_committed_time"]) + 0.5,
            int(checkpoint["cycle"]),
        )
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "cycle/time binding"):
            _history([self.passing], particle_time, particle_payload_overrides={5: payload})

    def test_ps_escape_requires_zero_initial_escape_source_closure_and_monotonic_ledgers(self) -> None:
        initial = _runtime_payload(self.passing)
        initial["ps_escape_checkpoints"][5]["ps_escape_ledger"][
            "ps_escaped_initial_cr_count_global"
        ] = 1.0
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "initial-CR escape"):
            _history([self.passing], initial)

        startup = _runtime_payload(self.passing)
        startup["ps_escape_checkpoints"][5]["ps_escape_ledger"][
            "ps_removed_excluded_early_cohort"
        ] = False
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "startup cohort removal"):
            _history([self.passing], startup)

        closure = _runtime_payload(self.passing)
        closure["ps_escape_checkpoints"][5]["active_injected_cr_count_global"] -= 1.0
        closure["ps_escape_checkpoints"][5]["active_injected_cr_mass_global"] -= 1.0
        closure["ps_escape_checkpoints"][5][
            "active_injected_cr_kinetic_energy_global"
        ] = _active_energy(999)
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError,
            "active \\+ startup_removed \\+ escaped_injected count",
        ):
            _history([self.passing], closure)

        nonmonotonic = _runtime_payload(self.passing)
        nonmonotonic["ps_escape_checkpoints"][5]["ps_escape_ledger"][
            "ps_injected_cr_energy_global"
        ] = 2.0
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "is not monotonic"):
            _history([self.passing], nonmonotonic)

    def test_one_cycle_self_declared_interval_endpoint_drift_and_gap_reject(self) -> None:
        one_cycle = _runtime_payload(self.passing)
        one_cycle["per_cycle_inventory"] = [
            {
                **one_cycle["per_cycle_inventory"][0],
                "start_time": 45.0,
                "end_time": 1200.0,
            }
        ]
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "too short"):
            _history([self.passing], one_cycle)
        endpoint = _runtime_payload(self.passing)
        endpoint["per_cycle_inventory"][-1]["end_time"] = 1199.5
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "terminal endpoint"):
            _history([self.passing], endpoint)
        gap = _runtime_payload(self.passing)
        gap["per_cycle_inventory"][1]["start_time"] += 0.1
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "time gap, overlap"):
            _history([self.passing], gap)

    def test_1155_invented_inline_cycle_records_and_arbitrary_source_are_non_admitting(
        self,
    ) -> None:
        invented = _runtime_payload(self.passing)
        invented["per_cycle_inventory"] = invented["per_cycle_inventory"][1:]
        self.assertEqual(len(invented["per_cycle_inventory"]), 1155)
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "binding keys drifted"):
            _history([self.passing], invented, bind_cycle_telemetry=False)

        arbitrary_source = _runtime_payload(self.passing)
        arbitrary_source["per_cycle_inventory"][500]["source_commit"] = "4" * 40
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "trusted source/escape identity drifted"
        ):
            _history([self.passing], arbitrary_source)

        fabricated_receipt = _runtime_payload(self.passing)
        fabricated_receipt["source_identity_receipt_sha256"] = "9" * 64
        for entry in fabricated_receipt["per_cycle_inventory"]:
            entry["source_identity_receipt_sha256"] = "9" * 64
        fabricated_receipt["escaped_slope_cutoff_evidence"][
            "source_identity_receipt_sha256"
        ] = "9" * 64
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "snapshot provenance disagrees"
        ):
            _history([self.passing], fabricated_receipt)

    def test_snapshot_raw_products_must_be_trusted_byte_decodable(self) -> None:
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "not a trusted Athena binary"
        ):
            _snapshot(ascii_raw_product="prtcl_jx")
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "lacks canonical prtcl_all execution metadata"
        ):
            _snapshot(ascii_raw_product="prtcl_all")

        tampered_after_reduction = _snapshot()
        raw_binding = tampered_after_reduction.record["snapshot_provenance"]["raw_products"][
            "mhd_w_bcc"
        ]
        (tampered_after_reduction.evidence_root / raw_binding["path"]).write_bytes(
            b"tampered after snapshot reduction\n"
        )
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "artifact byte count drifted"):
            _history([tampered_after_reduction], _runtime_payload(tampered_after_reduction))

    def test_history_recomputes_every_derived_snapshot_result_from_raw_bytes(self) -> None:
        attacks: list[app.ApplicabilitySnapshot] = []

        record = copy.deepcopy(dict(self.passing.record))
        maps = {name: np.array(values, copy=True) for name, values in self.passing.cell_maps.items()}
        maps["actual_leaf_dx1"][0, 0] += 1.0
        for values in maps.values():
            values.setflags(write=False)
        record["cell_map_contract"]["map_bindings"]["actual_leaf_dx1"]["sha256"] = (
            app._array_sha256(maps["actual_leaf_dx1"])
        )
        attacks.append(app.ApplicabilitySnapshot(record, maps, self.passing.evidence_root))

        record = copy.deepcopy(dict(self.passing.record))
        record["ion_scale_separation"]["precursor_magnetic_spectrum"][
            "resolved_restricted_power"
        ] += 1.0
        attacks.append(
            app.ApplicabilitySnapshot(record, self.passing.cell_maps, self.passing.evidence_root)
        )

        record = copy.deepcopy(dict(self.passing.record))
        record["particle_R_Lambda_exposure"]["sampling"] = "caller_fabricated_sampling"
        attacks.append(
            app.ApplicabilitySnapshot(record, self.passing.cell_maps, self.passing.evidence_root)
        )

        record = copy.deepcopy(dict(self.passing.record))
        record["gates"]["Q011-APP-NORM"]["threshold_provenance"] = "caller_fabricated_gate"
        attacks.append(
            app.ApplicabilitySnapshot(record, self.passing.cell_maps, self.passing.evidence_root)
        )

        for forged in attacks:
            with self.assertRaisesRegex(
                app.PhysicalApplicabilityError, "complete raw-byte recomputation"
            ):
                _history([forged], _runtime_payload(forged))

    def test_intermediate_symlink_path_escape_and_identity_receipts_reject(self) -> None:
        intermediate = _snapshot()
        raw = intermediate.evidence_root / "raw"
        raw.rename(intermediate.evidence_root / "raw-real")
        raw.symlink_to("raw-real", target_is_directory=True)
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "symlink"):
            _history([intermediate], _runtime_payload(intermediate))

        executable_alias = _snapshot()
        executable = executable_alias.evidence_root / "bound/athena"
        executable.rename(executable_alias.evidence_root / "bound/athena-real")
        executable.symlink_to("athena-real")
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "symlink"):
            _history([executable_alias], _runtime_payload(executable_alias))

        receipt_tamper = _snapshot()
        receipt = receipt_tamper.evidence_root / "bound/source-identity-receipt.json"
        receipt.write_bytes(b"caller-fabricated identity receipt\n")
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "receipt artifact .* drifted"):
            _history([receipt_tamper], _runtime_payload(receipt_tamper))

        escaped_path = copy.deepcopy(dict(self.passing.record))
        escaped_path["snapshot_provenance"]["raw_products"]["mhd_w_bcc"]["path"] = (
            "../escaped/mhd.bin"
        )
        forged = app.ApplicabilitySnapshot(
            escaped_path, self.passing.cell_maps, self.passing.evidence_root
        )
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "safe canonical relative"):
            _history([forged], _runtime_payload(forged))

    def test_checkpoint_must_be_first_committed_cycle_crossing_nominal_slot(self) -> None:
        payload = _runtime_payload(self.passing)
        checkpoint = payload["ps_escape_checkpoints"][0]
        first_crossing_index = next(
            index
            for index, entry in enumerate(payload["per_cycle_inventory"])
            if entry["end_time"] >= checkpoint["nominal_checkpoint_time"]
        )
        second = payload["per_cycle_inventory"][first_crossing_index + 1]
        checkpoint["observed_committed_time"] = second["end_time"]
        checkpoint["cycle"] = second["cycle"]
        checkpoint["previous_committed_cycle"] = second["previous_committed_cycle"]
        checkpoint["previous_committed_time"] = second["previous_committed_time"]
        checkpoint["ps_escape_ledger"] = _ps_escape_ledger(
            cycle=second["cycle"],
            observed_time=second["end_time"],
            active_count=1000,
            active_mass=1000.0,
        )
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError,
            "not the first committed cycle crossing",
        ):
            _history([self.passing], payload)

    def test_startup_coverage_requires_first_interval_crossing_and_no_post45_gap(
        self,
    ) -> None:
        passing = _history([self.passing], _runtime_payload(self.passing))
        coverage = passing["runtime_time_escape_evidence"]["cycle_coverage"]
        self.assertEqual(coverage["post_startup_removal_start_time"], 45.0)
        self.assertLess(coverage["first_interval_crossing_start_time"], 45.0)
        self.assertGreaterEqual(coverage["first_interval_crossing_end_time"], 45.0)
        self.assertLess(coverage["previous_committed_time_before_startup_crossing"], 45.0)

        bad_previous = _runtime_payload(self.passing)
        bad_previous["per_cycle_inventory"][0].update(
            {"previous_committed_time": 45.0, "start_time": 45.0, "end_time": 45.25}
        )
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "first committed interval crossing"
        ):
            _history([self.passing], bad_previous)

        omitted_crossing = _runtime_payload(self.passing)
        omitted_crossing["per_cycle_inventory"].pop(0)
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "first committed interval crossing"
        ):
            _history([self.passing], omitted_crossing)

        post_45_gap = _runtime_payload(self.passing)
        post_45_gap["per_cycle_inventory"][1]["start_time"] = 45.25
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "time gap"):
            _history([self.passing], post_45_gap)

    def test_snapshot_cycle_time_pair_must_be_an_admitted_telemetry_endpoint(self) -> None:
        shifted = _runtime_payload(self.passing)
        _shift_runtime_cycles(shifted, 2000)
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "admitted telemetry endpoint"):
            _history([self.passing], shifted)

    def test_any_inner_x1_escape_rejects_even_with_reason_code(self) -> None:
        payload = _runtime_payload(self.passing)
        payload["boundary_escape_ledger"]["nonperiodic_faces"]["inner_x1"] = _face_state(
            app.INNER_X1_ESCAPE_REASON, 1, 1.0, 1.0
        )
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "any inner_x1 escape"):
            _history([self.passing], payload)

        wrong_reason = _runtime_payload(self.passing)
        wrong_reason["boundary_escape_ledger"]["nonperiodic_faces"][
            "outer_x1"
        ]["reason_code"] = "arbitrary_escape_reason"
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "reason code drifted"):
            _history([self.passing], wrong_reason)

    def test_slope_cutoff_claim_requires_bound_complete_binwise_tail_evidence(self) -> None:
        unbound = _runtime_payload(self.passing)
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "binding keys drifted"):
            _history([self.passing], unbound, bind_slope_cutoff_evidence=False)

        unavailable = _runtime_payload(self.passing)
        _set_slope_cutoff_unavailable(unavailable)
        history = _history([self.passing], unavailable)
        applicability = history["runtime_time_escape_evidence"][
            "claim_specific_escape_applicability"
        ]
        self.assertFalse(applicability["high_energy_slope_or_cutoff_claim"]["pass"])
        self.assertTrue(applicability["Emax_claim"]["pass"])
        self.assertIn("high_energy_slope_or_cutoff_claim", history["claim_rejections"])

        biased = _runtime_payload(self.passing)
        ledger = biased["boundary_escape_ledger"]
        active_energy = _active_energy(999)
        ledger["nonperiodic_faces"]["outer_x1"] = _face_state(
            app.OUTER_X1_ESCAPE_REASON, 1, 1.0, 1.0
        )
        ledger["accumulated_escaped"] = _state_vector(1, 1.0, 1.0)
        ledger["terminal_active"] = _state_vector(999, 999.0, active_energy)
        last_entry = biased["per_cycle_inventory"][-1]
        last_entry["escaped_particle_specific_kinetic_energy_maximum"] = 1.0
        last_entry["escaped_particle_rg_maximum_over_Ly"] = 0.01
        last_entry["escaped_high_energy_tail_rg_maximum_over_Ly"] = 0.01
        last_entry["particle_rg_maximum_over_Ly"] = max(
            last_entry["particle_rg_maximum_over_Ly"], 0.01
        )
        last_entry["high_energy_tail_rg_maximum_over_Ly"] = max(
            last_entry["high_energy_tail_rg_maximum_over_Ly"], 0.01
        )
        terminal = biased["ps_escape_checkpoints"][-1]
        terminal["active_injected_cr_count_global"] = 999.0
        terminal["active_injected_cr_mass_global"] = 999.0
        terminal["active_injected_cr_kinetic_energy_global"] = active_energy
        terminal["escaped_injected_max_specific_kinetic_energy_global"] = 1.0
        terminal["escaped_injected_max_rg_over_Ly_global"] = 0.01
        terminal["ps_escape_ledger"] = _ps_escape_ledger(
            cycle=terminal["cycle"],
            observed_time=terminal["observed_committed_time"],
            active_count=999,
            active_mass=999.0,
            escaped_count=1,
            escaped_mass=1.0,
            escaped_energy=1.0,
        )
        _set_escape_events(
            biased,
            count=1,
            macro_weight=1.0,
            kinetic_energy=1.0,
            rg_over_ly=0.01,
            high_energy_tail_member=True,
        )
        evidence = biased["escaped_slope_cutoff_evidence"]
        (
            evidence["active_particle_count_by_bin"],
            evidence["active_macro_weight_by_bin"],
            evidence["active_kinetic_energy_by_bin"],
        ) = _active_energy_bins(999, evidence["energy_bin_edges"])
        evidence["escaped_particle_count_by_bin"] = [1, 0, 0, 0, 0]
        evidence["escaped_macro_weight_by_bin"] = [1.0, 0.0, 0.0, 0.0, 0.0]
        evidence["escaped_kinetic_energy_by_bin"] = [1.0, 0.0, 0.0, 0.0, 0.0]
        history = _history([self.passing], biased)
        slope = history["runtime_time_escape_evidence"]["slope_cutoff_escape_evidence"]
        self.assertFalse(slope["pass"])
        self.assertGreater(
            slope["escaped_fraction_by_bin"]["particle_count"][0],
            app.SLOPE_CUTOFF_BIN_ESCAPE_FRACTION_MAXIMUM,
        )
        self.assertFalse(
            history["runtime_time_escape_evidence"]["claim_specific_escape_applicability"][
                "high_energy_slope_or_cutoff_claim"
            ]["pass"]
        )

    def test_energy_bins_must_be_derived_from_terminal_PVTK_and_raw_escape_events(self) -> None:
        fabricated_active = _runtime_payload(self.passing)
        active_bins = fabricated_active["escaped_slope_cutoff_evidence"]
        active_bins["active_particle_count_by_bin"][0] -= 1
        active_bins["active_particle_count_by_bin"][1] += 1
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "not derived from bound raw particles/events"
        ):
            _history([self.passing], fabricated_active)

        fabricated_escape = _runtime_payload(self.passing)
        _configure_single_escape(fabricated_escape)
        escaped_bins = fabricated_escape["escaped_slope_cutoff_evidence"]
        escaped_bins["escaped_particle_count_by_bin"][0] = 0
        escaped_bins["escaped_particle_count_by_bin"][1] = 1
        escaped_bins["escaped_macro_weight_by_bin"][0] = 0.0
        escaped_bins["escaped_macro_weight_by_bin"][1] = 1.0
        escaped_bins["escaped_kinetic_energy_by_bin"][0] = 0.0
        escaped_bins["escaped_kinetic_energy_by_bin"][1] = 1.0
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "not derived from bound raw particles/events"
        ):
            _history([self.passing], fabricated_escape)

        unbound_escape = _runtime_payload(self.passing)
        _configure_single_escape(unbound_escape)
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "binding keys drifted"):
            _history([self.passing], unbound_escape, bind_escape_events=False)

        wrong_event_endpoint = _runtime_payload(self.passing)
        _configure_single_escape(wrong_event_endpoint)
        wrong_event_endpoint["escaped_particle_events"][0]["observed_committed_time"] -= 0.5
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "admitted telemetry endpoint"):
            _history([self.passing], wrong_event_endpoint)

    def test_escape_count_weight_energy_limits_and_momentum_residual_are_enforced(self) -> None:
        residual = _runtime_payload(self.passing)
        residual["boundary_escape_ledger"]["nonperiodic_faces"]["outer_x1"]["momentum"] = [
            1.0,
            0.0,
            0.0,
        ]
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "momentum residual"):
            _history([self.passing], residual)

        energetic_escape = _runtime_payload(self.passing)
        ledger = energetic_escape["boundary_escape_ledger"]
        ledger["nonperiodic_faces"]["outer_x1"] = _face_state(
            app.OUTER_X1_ESCAPE_REASON, 1, 1.0, 1000.0
        )
        ledger["accumulated_escaped"] = _state_vector(1, 1.0, 1000.0)
        ledger["terminal_active"] = _state_vector(999, 999.0, 1.0)
        last_entry = energetic_escape["per_cycle_inventory"][-1]
        last_entry["maximum_particle_specific_kinetic_energy"] = 1000.0
        last_entry["escaped_particle_specific_kinetic_energy_maximum"] = 1000.0
        last_entry["escaped_particle_rg_maximum_over_Ly"] = 0.01
        last_entry["particle_rg_maximum_over_Ly"] = max(
            last_entry["particle_rg_maximum_over_Ly"], 0.01
        )
        terminal = energetic_escape["ps_escape_checkpoints"][-1]
        terminal["active_injected_cr_count_global"] = 999.0
        terminal["active_injected_cr_mass_global"] = 999.0
        active_energy = _active_energy(999)
        ledger["terminal_active"] = _state_vector(999, 999.0, active_energy)
        terminal["active_injected_cr_kinetic_energy_global"] = active_energy
        terminal["escaped_injected_max_specific_kinetic_energy_global"] = 1000.0
        terminal["escaped_injected_max_rg_over_Ly_global"] = 0.01
        terminal["ps_escape_ledger"] = _ps_escape_ledger(
            cycle=terminal["cycle"],
            observed_time=terminal["observed_committed_time"],
            active_count=999,
            active_mass=999.0,
            escaped_count=1,
            escaped_mass=1.0,
            escaped_energy=1000.0,
        )
        _set_escape_events(
            energetic_escape,
            count=1,
            macro_weight=1.0,
            kinetic_energy=1000.0,
            rg_over_ly=0.01,
        )
        _set_slope_cutoff_unavailable(energetic_escape)
        history = _history([self.passing], energetic_escape)
        self.assertFalse(history["gates"]["Q011-APP-RG"]["pass"])
        self.assertGreater(
            history["gates"]["Q011-APP-RG"]["escape_fractions"]["kinetic_energy"],
            app.ESCAPED_KINETIC_ENERGY_FRACTION_MAXIMUM,
        )
        for claim, applicability in history["gates"]["Q011-APP-RG"][
            "claim_specific_escape_applicability"
        ].items():
            self.assertFalse(applicability["pass"], claim)
            self.assertIn(claim, history["claim_rejections"])
        self.assertIn("acceleration_efficiency_claim", history["claim_rejections"])

        checkpoint_maximum_drift = copy.deepcopy(energetic_escape)
        checkpoint_maximum_drift["ps_escape_checkpoints"][-1][
            "escaped_injected_max_specific_kinetic_energy_global"
        ] = 999.0
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "escaped maxima disagree"
        ):
            _history([self.passing], checkpoint_maximum_drift)

        missing_gyroradius = copy.deepcopy(energetic_escape)
        missing_gyroradius["per_cycle_inventory"][-1][
            "escaped_particle_rg_maximum_over_Ly"
        ] = 0.0
        missing_gyroradius["ps_escape_checkpoints"][-1][
            "escaped_injected_max_rg_over_Ly_global"
        ] = 0.0
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "raw-event maxima disagree"
        ):
            _history([self.passing], missing_gyroradius)

        count_mass_escape = _runtime_payload(self.passing)
        ledger = count_mass_escape["boundary_escape_ledger"]
        ledger["nonperiodic_faces"]["outer_x1"] = _face_state(
            app.OUTER_X1_ESCAPE_REASON, 2, 2.0, 1.0
        )
        ledger["accumulated_escaped"] = _state_vector(2, 2.0, 1.0)
        active_energy = _active_energy(998)
        ledger["terminal_active"] = _state_vector(998, 998.0, active_energy)
        last_entry = count_mass_escape["per_cycle_inventory"][-1]
        last_entry["escaped_particle_specific_kinetic_energy_maximum"] = 0.5
        last_entry["escaped_particle_rg_maximum_over_Ly"] = 0.01
        last_entry["particle_rg_maximum_over_Ly"] = max(
            last_entry["particle_rg_maximum_over_Ly"], 0.01
        )
        terminal = count_mass_escape["ps_escape_checkpoints"][-1]
        terminal["active_injected_cr_count_global"] = 998.0
        terminal["active_injected_cr_mass_global"] = 998.0
        terminal["active_injected_cr_kinetic_energy_global"] = active_energy
        terminal["escaped_injected_max_specific_kinetic_energy_global"] = 0.5
        terminal["escaped_injected_max_rg_over_Ly_global"] = 0.01
        terminal["ps_escape_ledger"] = _ps_escape_ledger(
            cycle=terminal["cycle"],
            observed_time=terminal["observed_committed_time"],
            active_count=998,
            active_mass=998.0,
            escaped_count=2,
            escaped_mass=2.0,
            escaped_energy=1.0,
        )
        _set_escape_events(
            count_mass_escape,
            count=2,
            macro_weight=2.0,
            kinetic_energy=1.0,
            rg_over_ly=0.01,
        )
        _set_slope_cutoff_unavailable(count_mass_escape)
        history = _history([self.passing], count_mass_escape)
        fractions = history["gates"]["Q011-APP-RG"]["escape_fractions"]
        self.assertGreater(
            fractions["particle_count"], app.ESCAPED_PARTICLE_COUNT_FRACTION_MAXIMUM
        )
        self.assertGreater(
            fractions["macro_weight"], app.ESCAPED_MACRO_WEIGHT_FRACTION_MAXIMUM
        )
        self.assertLess(
            fractions["kinetic_energy"], app.ESCAPED_KINETIC_ENERGY_FRACTION_MAXIMUM
        )
        for claim, applicability in history["gates"]["Q011-APP-RG"][
            "claim_specific_escape_applicability"
        ].items():
            self.assertFalse(applicability["pass"], claim)

    def test_escaped_particles_must_be_in_all_cycle_high_energy_and_gyroradius_extrema(self) -> None:
        omitted = _runtime_payload(self.passing)
        omitted["particle_exposure"][
            "escaped_particles_included_in_high_energy_and_gyroradius_extrema"
        ] = False
        history = _history([self.passing], omitted)
        self.assertFalse(history["gates"]["Q011-APP-RG"]["pass"])
        self.assertFalse(history["gates"]["Q011-APP-TIME"]["pass"])

        inconsistent = _runtime_payload(self.passing)
        inconsistent["per_cycle_inventory"][0]["escaped_particle_rg_maximum_over_Ly"] = 0.2
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "omits escaped particles"):
            _history([self.passing], inconsistent)

    def test_incomplete_escape_or_exposure_rejects_every_acceleration_claim(self) -> None:
        for section, field in (
            ("boundary_escape_ledger", "complete"),
            ("particle_exposure", "complete"),
        ):
            payload = _runtime_payload(self.passing)
            payload[section][field] = False
            history = _history([self.passing], payload)
            self.assertFalse(history["gates"]["Q011-APP-TIME"]["pass"])
            self.assertFalse(history["gates"]["Q011-APP-RG"]["pass"])
            for claim in (
                "Emax_claim",
                "high_energy_slope_or_cutoff_claim",
                "acceleration_rate_claim",
                "acceleration_efficiency_claim",
            ):
                self.assertIn(claim, history["claim_rejections"])

    def test_all_cycle_Lambda_or_DI_failure_rejects_acceleration_claims(self) -> None:
        for field, value, gate in (
            ("Lambda_maximum", 1.0, "Q011-APP-LAMBDA"),
            ("delta_B_rms_over_B0_minimum", 0.01, "Q011-APP-DI"),
        ):
            payload = _runtime_payload(self.passing)
            for entry in payload["per_cycle_inventory"]:
                entry[field] = value
            history = _history([self.passing], payload)
            self.assertFalse(history["gates"][gate]["pass"])
            for claim in (
                "Emax_claim",
                "high_energy_slope_or_cutoff_claim",
                "acceleration_rate_claim",
                "acceleration_efficiency_claim",
            ):
                self.assertIn(claim, history["claim_rejections"])

    def test_runtime_extrema_must_conservatively_dominate_retained_snapshots(self) -> None:
        payload = _runtime_payload(self.passing)
        for entry in payload["per_cycle_inventory"]:
            entry["S_delta_minimum_excluding_shock_transition"] += 1.0
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "do not conservatively dominate"):
            _history([self.passing], payload)

    def test_readiness_record_binds_only_successor_files_and_no_fabricated_telemetry(self) -> None:
        readiness = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertEqual(readiness["successor_id"], app.SUCCESSOR_ID)
        self.assertEqual(readiness["production_evidence_status"], "absent_not_fabricated")
        self.assertFalse(readiness["runtime_time_escape_telemetry_available"])
        for binding in readiness["source_bindings"].values():
            path = REPO_ROOT / binding["path"]
            self.assertEqual(binding["sha256"], _sha256(path))
        self.assertEqual(
            readiness["source_bindings"]["design"]["path"], str(DESIGN.relative_to(REPO_ROOT))
        )


if __name__ == "__main__":
    unittest.main()
