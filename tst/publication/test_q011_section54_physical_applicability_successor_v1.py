#!/usr/bin/env python3
"""Focused adversarial tests for the hardened Q011 applicability successor."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
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
SOURCE_COMMIT = "4" * 40
ATTEMPT_ID = "q011-applicability-adversarial-fixture"
RAW_PATHS = {
    product: f"raw/{product}.00500.bin"
    for product in app.REQUIRED_RAW_PRODUCTS
}
RAW_PATHS["prtcl_all"] = "raw/prtcl_all.00500.part.vtk"


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
        (0.0, NX1 * DX1, 0.0, NX2 * DX2, 0.0, 1.0),
        fields,
    )
    return output_primitives.AthenaBinaryDataset(
        source=RAW_PATHS["mhd_w_bcc"],
        time=TIME,
        cycle=500,
        location_size=8,
        variable_size=8,
        variable_names=tuple(reversed(science.MHD_PRIMITIVE_FIELDS)),
        input_parameters={},
        root_grid_shape=(NX1, NX2, 1),
        meshblock_shape=(NX1, NX2, 1),
        nghost=0,
        domain_bounds=(0.0, NX1 * DX1, 0.0, NX2 * DX2, 0.0, 1.0),
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
        input_parameters={},
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


def _snapshot_evidence(
    root: Path,
    mhd: output_primitives.AthenaBinaryDataset,
    currents: dict[str, output_primitives.AthenaBinaryDataset],
    particles: dict[str, object],
    *,
    normalization: object | None = None,
    tamper_deck_after_binding: bool = False,
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
    raw = {
        product: _write_artifact(root, product, RAW_PATHS[product], product.encode("ascii"))
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
        "source_archive_sha256": bindings["source_archive"]["sha256"],
        "executable_sha256": bindings["executable"]["sha256"],
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
) -> app.ApplicabilitySnapshot:
    with tempfile.TemporaryDirectory() as temporary:
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
        )
        if tamper_decoded_particle_after_manifest:
            particles["velocity"][0, 0] += 1.0
        return app.reduce_physical_applicability_snapshot(
            mhd,
            currents,
            evidence_root=root,
            normalization_evidence=normalization_evidence,
            snapshot_provenance=provenance,
            particle_source=RAW_PATHS["prtcl_all"],
            nominal_slot_time=TIME,
            observed_committed_time=TIME,
            **particles,
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


def _active_energy(count: int, velocity: float = 10.0) -> float:
    return count * _specific_energy_bounds(velocity)[1]


def _particle_checkpoint_vtk(
    count: int, observed_time: float, cycle: int, *, velocity: float = 10.0
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
    velocities[:, 0] = velocity
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
    inventory = [
        {
            "cycle": 1000 + index,
            "start_time": app.STARTUP_REMOVAL_TIME + index,
            "end_time": app.STARTUP_REMOVAL_TIME + index + 1.0,
            "telemetry_record_sha256": _sha256_bytes(f"telemetry-{index}".encode("ascii")),
            **metrics,
        }
        for index in range(int(app.EXPECTED_TERMINAL_TIME - app.STARTUP_REMOVAL_TIME))
    ]
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
    return {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.RUNTIME_RECORD_TYPE,
        "successor_id": app.SUCCESSOR_ID,
        "qualification_effect": app.QUALIFICATION_EFFECT,
        "authorization": dict(app.AUTHORIZATION),
        "attempt_id": record["snapshot_provenance"]["attempt_id"],
        "source_commit": record["snapshot_provenance"]["source_commit"],
        "executable_sha256": record["snapshot_provenance"]["executable_sha256"],
        "runtime_normalization_sha256": record["snapshot_provenance"][
            "runtime_normalization_sha256"
        ],
        "ps_escape_accounting_source_commit": app.PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
        "per_cycle_inventory": inventory,
        "ps_escape_checkpoints": checkpoints,
        "claim_escape_accounting": dict(app._CLAIM_ESCAPE_ACCOUNTING),
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
                "ix1": _state_vector(),
                "ox1": _state_vector(),
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


def _history(
    snapshots: list[app.ApplicabilitySnapshot],
    payload: dict[str, object],
    *,
    restart_ledger_overrides: dict[int, dict[str, object]] | None = None,
    particle_count_overrides: dict[int, int] | None = None,
    particle_payload_overrides: dict[int, bytes] | None = None,
) -> dict[str, object]:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        bound_payload = copy.deepcopy(payload)
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
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "decoded prtcl_all payload"):
            _snapshot(tamper_decoded_particle_after_manifest=True)
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "unexpected keyword"):
            app.reduce_physical_applicability_snapshot(  # type: ignore[call-arg]
                None, {}, normalization=dict(app.EXACT_NORMALIZATION)
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
            stripped = app.ApplicabilitySnapshot(record=record, cell_maps=self.passing.cell_maps)
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
        self.assertEqual(coverage["covered_cycle_count"], 1155)
        self.assertEqual(coverage["post_startup_removal_start_time"], 45.0)
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
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "exact actual endpoints"):
            _history([self.passing], endpoint)
        gap = _runtime_payload(self.passing)
        gap["per_cycle_inventory"][1]["start_time"] += 0.1
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "time gap or overlap"):
            _history([self.passing], gap)

    def test_escape_count_weight_energy_limits_and_momentum_residual_are_enforced(self) -> None:
        residual = _runtime_payload(self.passing)
        residual["boundary_escape_ledger"]["nonperiodic_faces"]["ix1"]["momentum"] = [
            1.0,
            0.0,
            0.0,
        ]
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "momentum residual"):
            _history([self.passing], residual)

        energetic_escape = _runtime_payload(self.passing)
        ledger = energetic_escape["boundary_escape_ledger"]
        ledger["nonperiodic_faces"]["ix1"] = _state_vector(1, 1.0, 1000.0)
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
            app.PhysicalApplicabilityError, "positive escaped census lacks"
        ):
            _history([self.passing], missing_gyroradius)

        count_mass_escape = _runtime_payload(self.passing)
        ledger = count_mass_escape["boundary_escape_ledger"]
        ledger["nonperiodic_faces"]["ix1"] = _state_vector(2, 2.0, 1.0)
        ledger["accumulated_escaped"] = _state_vector(2, 2.0, 1.0)
        active_energy = _active_energy(998)
        ledger["terminal_active"] = _state_vector(998, 998.0, active_energy)
        last_entry = count_mass_escape["per_cycle_inventory"][-1]
        last_entry["escaped_particle_specific_kinetic_energy_maximum"] = 1.0
        last_entry["escaped_particle_rg_maximum_over_Ly"] = 0.01
        last_entry["particle_rg_maximum_over_Ly"] = max(
            last_entry["particle_rg_maximum_over_Ly"], 0.01
        )
        terminal = count_mass_escape["ps_escape_checkpoints"][-1]
        terminal["active_injected_cr_count_global"] = 998.0
        terminal["active_injected_cr_mass_global"] = 998.0
        terminal["active_injected_cr_kinetic_energy_global"] = active_energy
        terminal["escaped_injected_max_specific_kinetic_energy_global"] = 1.0
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
