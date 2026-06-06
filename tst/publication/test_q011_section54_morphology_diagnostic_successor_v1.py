#!/usr/bin/env python3
"""Focused adversarial tests for the Q011 t=500 morphology successor."""

from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

try:
    from tst.publication import analyze_q011_section54_outputs as output_primitives
    from tst.publication import q011_section54_morphology_diagnostic_successor_v1 as morphology
    from tst.publication import q011_section54_production_science_successor_v1 as science
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_morphology_diagnostic_successor_v1 as morphology
    import q011_section54_production_science_successor_v1 as science


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_morphology_diagnostic_successor_v1_2026-06-06.json"
)
ROOT_SHAPE = (700, 4, 1)
BLOCK_SHAPE = (350, 4, 1)
DOMAIN = (0.0, 16800.0, 0.0, 96.0, 0.0, 1.0)
MIXED_LEAVES = (
    ((0, 0, 0), 1),
    ((0, 1, 0), 1),
    ((1, 0, 0), 1),
    ((1, 1, 0), 1),
    ((1, 0, 0), 0),
)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json(payload: bytes) -> dict[str, object]:
    return json.loads(payload.decode("utf-8"))


def _geometry(
    logical_location: tuple[int, int, int], level: int
) -> tuple[float, float, float, float, float, float]:
    result = []
    for axis, logical_index in enumerate(logical_location):
        root_blocks = ROOT_SHAPE[axis] // BLOCK_SHAPE[axis]
        logical_extent = root_blocks * (2**level if ROOT_SHAPE[axis] > 1 else 1)
        lower = DOMAIN[2 * axis]
        upper = DOMAIN[2 * axis + 1]
        width = (upper - lower) / logical_extent
        result.extend((lower + logical_index * width, lower + (logical_index + 1) * width))
    return tuple(result)


def _spatial_fields(
    geometry: tuple[float, float, float, float, float, float]
) -> dict[str, np.ndarray]:
    x1 = np.linspace(
        geometry[0], geometry[1], BLOCK_SHAPE[0], endpoint=False, dtype=np.float64
    )
    x1 += 0.5 * (geometry[1] - geometry[0]) / BLOCK_SHAPE[0]
    x2 = np.linspace(
        geometry[2], geometry[3], BLOCK_SHAPE[1], endpoint=False, dtype=np.float64
    )
    x2 += 0.5 * (geometry[3] - geometry[2]) / BLOCK_SHAPE[1]
    front_x = 4998.0
    dens_x = np.where(x1 < front_x, 4.0, 1.0)
    if np.any(x1 == front_x):
        dens_x[x1 == front_x] = 3.8
    upstream = x1 > front_x
    bcc1 = np.where(
        upstream[None, :],
        np.where(x2[:, None] < 48.0, 2.0, 4.0),
        1.0,
    )
    velx = np.broadcast_to(
        np.where(upstream[None, :], -30.0, 0.0),
        (BLOCK_SHAPE[1], BLOCK_SHAPE[0]),
    )
    dens = np.broadcast_to(dens_x[None, :], (BLOCK_SHAPE[1], BLOCK_SHAPE[0]))
    zeros = np.zeros_like(dens)
    return {
        "dens": dens[None, :, :],
        "velx": velx[None, :, :],
        "vely": zeros[None, :, :],
        "velz": zeros[None, :, :],
        "eint": np.full((1, BLOCK_SHAPE[1], BLOCK_SHAPE[0]), 1.5),
        "bcc1": bcc1[None, :, :],
        "bcc2": zeros[None, :, :],
        "bcc3": zeros[None, :, :],
    }


def _block(
    logical_location: tuple[int, int, int],
    level: int,
    fields: dict[str, np.ndarray],
) -> output_primitives.AthenaBinaryBlock:
    return output_primitives.AthenaBinaryBlock(
        index_bounds=(
            0,
            BLOCK_SHAPE[0] - 1,
            0,
            BLOCK_SHAPE[1] - 1,
            0,
            BLOCK_SHAPE[2] - 1,
        ),
        logical_location=logical_location,
        level=level,
        geometry=_geometry(logical_location, level),
        fields=fields,
    )


def _mhd_dataset(source: str) -> output_primitives.AthenaBinaryDataset:
    blocks = []
    for logical_location, level in MIXED_LEAVES:
        geometry = _geometry(logical_location, level)
        blocks.append(_block(logical_location, level, _spatial_fields(geometry)))
    return output_primitives.AthenaBinaryDataset(
        source=source,
        time=500.0,
        cycle=500,
        location_size=8,
        variable_size=8,
        variable_names=tuple(reversed(science.MHD_PRIMITIVE_FIELDS)),
        input_parameters={},
        root_grid_shape=ROOT_SHAPE,
        meshblock_shape=BLOCK_SHAPE,
        nghost=0,
        domain_bounds=DOMAIN,
        blocks=tuple(blocks),
    )


def _current_dataset(
    product: str, field: str, source: str
) -> output_primitives.AthenaBinaryDataset:
    value = {
        "prtcl_rho": 0.1,
        "prtcl_jx": 1.0,
        "prtcl_jy": 2.0,
        "prtcl_jz": 3.0,
        "mhd_j2": 14.0,
    }[product]
    blocks = tuple(
        _block(
            logical_location,
            level,
            {field: np.full((1, BLOCK_SHAPE[1], BLOCK_SHAPE[0]), value)},
        )
        for logical_location, level in MIXED_LEAVES
    )
    return output_primitives.AthenaBinaryDataset(
        source=source,
        time=500.0,
        cycle=500,
        location_size=8,
        variable_size=8,
        variable_names=(field,),
        input_parameters={},
        root_grid_shape=ROOT_SHAPE,
        meshblock_shape=BLOCK_SHAPE,
        nghost=0,
        domain_bounds=DOMAIN,
        blocks=blocks,
    )


def _velocity_from_chi(values: np.ndarray) -> np.ndarray:
    momentum = science.UPSTREAM_SPEED_U0 * np.sqrt(values)
    velocity = momentum / np.sqrt(
        1.0 + (momentum / science.PARTICLE_LIGHT_SPEED) ** 2
    )
    return np.column_stack((velocity, np.zeros_like(velocity), np.zeros_like(velocity))).astype(
        np.float32
    ).astype(np.float64)


def _particle_payload() -> dict[str, np.ndarray]:
    downstream_count = science.DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES
    upstream_count = 100
    count = downstream_count + upstream_count
    x1 = np.concatenate(
        (
            np.full(downstream_count, 4500.0),
            np.full(upstream_count, 5500.0),
        )
    )
    x2 = (6.0 + 12.0 * (np.arange(count) % 8)).astype(np.float64)
    points = np.column_stack((x1, x2, np.full(count, 0.5))).astype(np.float32).astype(
        np.float64
    )
    return {
        "points": points,
        "cr_source": np.ones(count, dtype=np.int64),
        "birth_time": np.full(count, 45.0),
        "velocity": _velocity_from_chi(np.linspace(1.0, 16.0, count)),
        "macro_weight": np.linspace(0.5, 1.5, count),
    }


def _raw_bindings() -> dict[str, dict[str, object]]:
    result = {}
    for product in morphology.EXPECTED_RAW_PRODUCTS:
        suffix = ".part.vtk" if product == "prtcl_all" else ".bin"
        directory = "pvtk" if product == "prtcl_all" else "bin"
        path = f"raw/{directory}/q011.{product}.00500{suffix}"
        result[product] = {
            "product": product,
            "path": path,
            "sha256": _sha256(product.encode("ascii")),
            "nominal_slot_time": 500.0,
            "observed_committed_time": 500.0,
            "cycle": 500,
        }
    return result


def _particle_metadata(bindings: dict[str, dict[str, object]]) -> dict[str, object]:
    return {
        "source": bindings["prtcl_all"]["path"],
        "nominal_slot_time": 500.0,
        "observed_committed_time": 500.0,
        "cycle": 500,
        "scalar_fields": list(morphology.EXPECTED_PVTK_SCALARS),
        "vector_fields": list(morphology.EXPECTED_PVTK_VECTORS),
        "storage_layout": "single_mpi_io_pvtk_file_with_nranks_header",
    }


def _inputs() -> dict[str, object]:
    bindings = _raw_bindings()
    mhd = _mhd_dataset(str(bindings["mhd_w_bcc"]["path"]))
    currents = {
        product: _current_dataset(product, field, str(bindings[product]["path"]))
        for product, field in science.CURRENT_PRODUCT_FIELDS.items()
    }
    return {
        "attempt_id": "q011-morphology-fixture-amr-seed-23050101",
        "mhd_dataset": mhd,
        "current_datasets": currents,
        "raw_artifact_bindings": bindings,
        "particle_snapshot_metadata": _particle_metadata(bindings),
        "observed_committed_time": 500.0,
        "cycle": 500,
        **_particle_payload(),
    }


class Q011Section54MorphologyDiagnosticSuccessorV1Tests(unittest.TestCase):
    def test_builds_deterministic_bound_profiles_and_refinement_overlay(self) -> None:
        first = morphology.build_morphology_artifacts(**_inputs())
        second = morphology.build_morphology_artifacts(**_inputs())
        self.assertEqual(first, second)
        self.assertEqual(tuple(first), morphology.ALL_MEMBERS)

        manifest = morphology.validate_morphology_artifacts(first)
        self.assertEqual(manifest["authorization"], dict(morphology.AUTHORIZATION))
        self.assertEqual(
            manifest["qualification_effect"], morphology.QUALIFICATION_EFFECT
        )

        particle = _json(first[morphology.PARTICLE_PROFILE_MEMBER])
        self.assertAlmostEqual(
            particle["normalization_closure"]["energy_fraction_sum"], 1.0
        )
        self.assertAlmostEqual(
            particle["normalization_closure"]["probability_density_integral"], 1.0
        )
        self.assertEqual(
            particle["normalization_closure"]["active_particle_count"],
            science.DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES + 100,
        )

        current = _json(first[morphology.CURRENT_PROFILE_MEMBER])
        self.assertIn("deposited_J_CR_over_c", current["source_representation"]["prtcl_j_components"])
        self.assertEqual(
            current["frame_transform"]["formula"],
            "(J_CR/c)_gas = (J_CR/c)_lab - (rho_CR/c) * u_gas",
        )
        self.assertEqual(
            len(current["profiles"]["deposited_j_cr_over_c_lab_x"]),
            int((DOMAIN[1] - DOMAIN[0]) / 12.0),
        )

        magnetic = _json(first[morphology.MAGNETIC_PROFILE_MEMBER])
        self.assertIn("gas_density_context", magnetic["profiles"])
        self.assertIn("bmag_over_b0", magnetic["profiles"])

        refinement = _json(first[morphology.REFINEMENT_PROFILE_MEMBER])
        np.testing.assert_allclose(refinement["area_fraction_sum_by_x1"], 1.0)
        self.assertEqual(
            set(refinement["profiles"]),
            {"source_level_0_area_fraction", "source_level_1_area_fraction"},
        )

        overlay = _json(first[morphology.REFINEMENT_OVERLAY_MEMBER])
        self.assertEqual(len(overlay["leaf_block_rectangles"]), len(MIXED_LEAVES))
        self.assertEqual(
            overlay["leaf_block_count_by_source_level"],
            [
                {"leaf_block_count": 1, "source_level": 0},
                {"leaf_block_count": 4, "source_level": 1},
            ],
        )

    def test_raw_and_decoded_provenance_relabeling_fails_closed(self) -> None:
        source_bindings = morphology.source_bindings()
        source_bindings.pop(next(iter(source_bindings)))
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "source binding inventory drifted"
        ):
            morphology._validate_source_binding_map(source_bindings)

        values = _inputs()
        values["raw_artifact_bindings"].pop("mhd_j2")
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "inventory drifted"
        ):
            morphology.build_morphology_artifacts(**values)

        values = _inputs()
        values["raw_artifact_bindings"]["prtcl_jx"]["path"] = (
            "raw/bin/q011.prtcl_jy.00500.bin"
        )
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "product identity"
        ):
            morphology.build_morphology_artifacts(**values)

        values = _inputs()
        values["current_datasets"]["prtcl_jx"] = replace(
            values["current_datasets"]["prtcl_jx"], source="raw/bin/relabel.prtcl_jx.00500.bin"
        )
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "decoded source"
        ):
            morphology.build_morphology_artifacts(**values)

    def test_particle_metadata_time_and_inventory_drift_fail_closed(self) -> None:
        values = _inputs()
        values["particle_snapshot_metadata"]["scalar_fields"].pop()
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "scalar inventory drifted"
        ):
            morphology.build_morphology_artifacts(**values)

        values = _inputs()
        values["particle_snapshot_metadata"]["cycle"] = 499
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "metadata cycle drifted"
        ):
            morphology.build_morphology_artifacts(**values)

        values = _inputs()
        values["raw_artifact_bindings"]["prtcl_all"]["nominal_slot_time"] = 400.0
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "must bind nominal t=500"
        ):
            morphology.build_morphology_artifacts(**values)

    def test_particle_domain_float32_and_energy_requirements_fail_closed(self) -> None:
        values = _inputs()
        values["points"][0, 0] = 4500.1
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "exactly representable decoded float32"
        ):
            morphology.build_morphology_artifacts(**values)

        values = _inputs()
        values["points"][0, 1] = 120.0
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "escaped retained x2 domain"
        ):
            morphology.build_morphology_artifacts(**values)

        values = _inputs()
        values["velocity"][:] = 0.0
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "kinetic energy must be positive"
        ):
            morphology.build_morphology_artifacts(**values)

    def test_misaligned_fixed_grid_and_artifact_tamper_fail_closed(self) -> None:
        values = _inputs()
        state = science.compose_full_mhd_state(values["mhd_dataset"])
        misaligned = replace(
            state,
            x1_faces=np.linspace(
                float(state.x1_faces[0]),
                float(state.x1_faces[-1]) + 1.0,
                state.x1_faces.size,
            ),
        )
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "fixed 12 c/omega_pi grid"
        ):
            morphology._restriction_shape(misaligned)

        artifacts = morphology.build_morphology_artifacts(**values)
        provenance_tampered = dict(artifacts)
        current = _json(provenance_tampered[morphology.CURRENT_PROFILE_MEMBER])
        current["provenance"]["source_products"].pop()
        provenance_tampered[morphology.CURRENT_PROFILE_MEMBER] = (
            morphology._canonical_json_bytes(current)
        )
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "source-product provenance"
        ):
            morphology.validate_morphology_artifacts(provenance_tampered)

        tampered = dict(artifacts)
        particle = _json(tampered[morphology.PARTICLE_PROFILE_MEMBER])
        particle["profiles"]["active_cr_kinetic_energy_fraction"][0] = 0.25
        tampered[morphology.PARTICLE_PROFILE_MEMBER] = morphology._canonical_json_bytes(
            particle
        )
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "inventory checksum"
        ):
            morphology.validate_morphology_artifacts(tampered)

        extra = {**artifacts, "unexpected.json": b"{}\n"}
        with self.assertRaisesRegex(
            morphology.MorphologyDiagnosticError, "member inventory drifted"
        ):
            morphology.validate_morphology_artifacts(extra)

    def test_readiness_binds_exact_sources_and_refuses_authority(self) -> None:
        contract = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertNotIn("REPLACE_", READINESS.read_text(encoding="utf-8"))
        self.assertEqual(contract["schema_version"], 1)
        self.assertEqual(contract["successor_id"], morphology.SUCCESSOR_ID)
        self.assertEqual(contract["authorization"], dict(morphology.AUTHORIZATION))
        self.assertEqual(
            contract["qualification_effect"], morphology.QUALIFICATION_EFFECT
        )
        for binding in contract["source_bindings"].values():
            path = REPO_ROOT / binding["path"]
            self.assertTrue(path.is_file())
            self.assertEqual(binding["sha256"], _sha256(path.read_bytes()))
        self.assertEqual(
            contract["artifact_contract"]["members"],
            list(morphology.ALL_MEMBERS),
        )
        self.assertEqual(
            contract["profile_contract"]["fixed_cell_size_c_over_omega_pi"],
            morphology.FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI,
        )


if __name__ == "__main__":
    unittest.main()
