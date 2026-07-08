#!/usr/bin/env python3
"""Focused tests for the corrected Bell deposited-current raw-output oracle."""

from __future__ import annotations

import copy
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import tempfile
import unittest

import numpy as np

from tst.publication import q043_bell_current_volume_aware_deposited_current_oracle as oracle


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q043_bell_current_volume_aware_deposited_current_oracle_source_local_2026-06-06.json"
)
PGEN_SOURCE = REPO_ROOT / "src/pgen/tests/q043_bell_current_volume_aware.cpp"
SOURCE_CONTRACT_HARNESS = (
    REPO_ROOT / "tst/publication/q043_bell_current_volume_aware_host_harness.cpp"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _case(case_id: str) -> dict[str, object]:
    return next(case for case in oracle.expected_cases() if case["case_id"] == case_id)


def _with_species0_overrides(text: str, overrides: dict[str, str]) -> str:
    additions = "".join(f"{name} = {value}\n" for name, value in overrides.items())
    marker = "\n<problem>\n"
    if text.count(marker) != 1:
        raise AssertionError("rendered deck species0 insertion point drifted")
    return text.replace(marker, f"{additions}{marker}", 1)


def _block_payload(
    case: dict[str, object],
    field: str,
    logical_location: tuple[int, int, int],
    *,
    values: np.ndarray,
) -> bytes:
    nx1, nx2, nx3 = (int(value) for value in case["meshblock_nx"])
    bounds = tuple(tuple(float(value) for value in item) for item in case["bounds"])
    meshblock_grid = tuple(int(value) for value in case["meshblock_grid"])
    geometry = tuple(
        bound
        for axis, (lower, upper) in enumerate(bounds)
        for bound in (
            lower
            + logical_location[axis] * (upper - lower) / meshblock_grid[axis],
            lower
            + (logical_location[axis] + 1)
            * (upper - lower)
            / meshblock_grid[axis],
        )
    )
    expected_shape = (nx3, nx2, nx1)
    if values.shape != expected_shape:
        raise AssertionError(f"{field}: expected shape {expected_shape}, got {values.shape}")
    index_and_logical = (
        0,
        nx1 - 1,
        0,
        nx2 - 1,
        0,
        nx3 - 1,
        *logical_location,
        0,
    )
    return (
        struct.pack("<10i", *index_and_logical)
        + struct.pack("<6d", *geometry)
        + np.asarray(values, dtype="<f4").tobytes()
    )


def _runtime_parameter_header(
    case: dict[str, object],
    overrides: dict[tuple[str, str], str] | None = None,
) -> str:
    overrides = overrides or {}
    parameters = oracle.expected_runtime_parameters(case)
    for (block, key), value in overrides.items():
        if block not in parameters:
            raise AssertionError(f"unknown runtime parameter block: {block}")
        if key not in parameters[block] and key not in oracle.RUNTIME_OUTPUT_BOOKKEEPING_KEYS:
            raise AssertionError(f"unknown runtime parameter override: {block}/{key}")
        parameters[block][key] = value
    return "\n".join(
        line
        for block, values in parameters.items()
        for line in (f"<{block}>", *(f"{key} = {value}" for key, value in values.items()))
    )


def _binary_payload(
    case: dict[str, object],
    field: str,
    blocks: list[bytes],
    *,
    parameter_overrides: dict[tuple[str, str], str] | None = None,
) -> bytes:
    parameter_header = _runtime_parameter_header(case, parameter_overrides).encode("utf-8")
    return (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        b"  time=1.0e-6\n"
        b"  cycle=1\n"
        b"  size of location=8\n"
        b"  size of variable=4\n"
        b"  number of variables=1\n"
        + f"  variables:  {field}  \n".encode("ascii")
        + f"  header offset={len(parameter_header)}\n".encode("ascii")
        + parameter_header
        + b"".join(blocks)
    )


def _field_value(case: dict[str, object], field: str) -> float:
    if field == "prtcl_rho":
        return oracle.EXPECTED_RHO
    basis = oracle._mode_basis(int(case["dimension"]))
    return oracle.EXPECTED_J_OVER_C * basis[oracle.FIELDS.index(field) - 1]


def _real_cycle_one_output_bookkeeping(field: str) -> dict[tuple[str, str], str]:
    field_index = oracle.FIELDS.index(field) + 1
    bookkeeping: dict[tuple[str, str], str] = {}
    for output_index in range(1, len(oracle.FIELDS) + 1):
        bookkeeping[(f"output{output_index}", "file_number")] = str(
            1 if output_index <= field_index else 2
        )
        bookkeeping[(f"output{output_index}", "last_time")] = "0"
    return bookkeeping


def _write_raw_case(
    root: Path,
    case: dict[str, object],
    *,
    mutations: dict[tuple[str, int], tuple[int, float]] | None = None,
    parameter_overrides: dict[str, dict[tuple[str, str], str]] | None = None,
) -> dict[str, tuple[Path, ...]]:
    mutations = mutations or {}
    parameter_overrides = parameter_overrides or {}
    paths: dict[str, tuple[Path, ...]] = {}
    shard_count = int(case["mpi_ranks"])
    shape = tuple(reversed(tuple(int(value) for value in case["meshblock_nx"])))
    blocks_x1, blocks_x2, blocks_x3 = (
        int(value) for value in case["meshblock_grid"]
    )
    logical_locations = tuple(
        (
            shard % blocks_x1,
            (shard // blocks_x1) % blocks_x2,
            shard // (blocks_x1 * blocks_x2),
        )
        for shard in range(blocks_x1 * blocks_x2 * blocks_x3)
    )
    if len(logical_locations) != shard_count:
        raise AssertionError("one-MeshBlock-per-rank fixture contract drifted")
    for field in oracle.FIELDS:
        field_paths = []
        runtime_overrides = _real_cycle_one_output_bookkeeping(field)
        runtime_overrides.update(parameter_overrides.get(field, {}))
        for shard in range(shard_count):
            values = np.full(shape, _field_value(case, field), dtype=np.float64)
            mutation = mutations.get((field, shard))
            if mutation is not None:
                flat_index, delta = mutation
                values.reshape(-1)[flat_index] += delta
            payload = _binary_payload(
                case,
                field,
                [_block_payload(case, field, logical_locations[shard], values=values)],
                parameter_overrides=runtime_overrides,
            )
            path = root / str(case["case_id"]) / field / f"rank{shard}.bin"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            field_paths.append(path)
        paths[field] = tuple(field_paths)
    return paths


def _write_raw_matrix(root: Path) -> dict[str, dict[str, tuple[Path, ...]]]:
    return {
        str(case["case_id"]): _write_raw_case(root, case)
        for case in oracle.expected_cases()
    }


class Q043BellCurrentVolumeAwareDepositedCurrentOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary_directory = tempfile.TemporaryDirectory()
        cls.source_contract_binary = (
            Path(cls._temporary_directory.name) / "q043-deck-source-contract"
        )
        subprocess.run(
            [
                os.environ.get("CXX", "c++"),
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(SOURCE_CONTRACT_HARNESS),
                "-o",
                str(cls.source_contract_binary),
            ],
            cwd=REPO_ROOT,
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary_directory.cleanup()

    def test_exact_deck_matrix_derives_qscale_from_root_volume_and_ppc(self) -> None:
        cases = oracle.expected_cases()
        self.assertEqual(len(cases), 132)
        self.assertEqual({case["dimension"] for case in cases}, {1, 2, 3})
        self.assertEqual({case["resolution"] for case in cases}, {"coarse", "fine"})
        self.assertEqual({case["ppc"] for case in cases}, {1, 4})
        self.assertEqual(
            {case["decomposition"] for case in cases},
            {"single", "split_x1", "split_x2", "split_x1x2", "split_x3", "split_xyz"},
        )
        self.assertEqual(
            {case["artificial_c_over_v_cr"] for case in cases},
            {100, 1000, 10000},
        )
        self.assertEqual({case["species_mass"] for case in cases}, {2.0})
        self.assertEqual({case["species_charge"] for case in cases}, {4.0 * math.pi * 1.0e-6})
        self.assertNotIn(800000.0, {case["deposit_qscale"] for case in cases})

        grouped: dict[tuple[int, str, int], list[dict[str, object]]] = {}
        for case in cases:
            grouped.setdefault(
                (int(case["dimension"]), str(case["resolution"]), int(case["ppc"])),
                [],
            ).append(case)
            expected = (
                oracle.EXPECTED_J_OVER_C
                * float(case["root_cell_volume"])
                / (int(case["ppc"]) * oracle.SPECIES_CHARGE * oracle.STREAM_SPEED)
            )
            self.assertAlmostEqual(float(case["deposit_qscale"]), expected)
            self.assertAlmostEqual(
                oracle.configured_volume_mean_j_over_c(case),
                oracle.EXPECTED_J_OVER_C,
            )
        for (dimension, _, _), pair in grouped.items():
            self.assertEqual(
                len(pair),
                len(oracle.DECOMPOSITIONS_BY_DIMENSION[dimension])
                * len(oracle.ARTIFICIAL_C_OVER_V_CR_VALUES),
            )
            self.assertEqual({item["deposit_qscale"] for item in pair}, {pair[0]["deposit_qscale"]})
        for case in cases:
            self.assertEqual(
                int(case["mpi_ranks"]),
                math.prod(int(value) for value in case["meshblock_grid"]),
            )
            self.assertEqual(
                tuple(int(value) for value in case["meshblock_nx"]),
                tuple(
                    int(count) // int(blocks)
                    for count, blocks in zip(case["global_nx"], case["meshblock_grid"])
                ),
            )
            self.assertTrue(
                all(int(count) == 1 or int(count) >= 4 for count in case["meshblock_nx"])
            )

    def test_artificial_c_and_species_mass_are_absent_from_deposited_current_formula(
        self,
    ) -> None:
        signature = inspect.signature(oracle.required_deposit_qscale)
        self.assertNotIn("artificial_light_speed", signature.parameters)
        self.assertNotIn("species_mass", signature.parameters)
        self.assertIn("species_charge", signature.parameters)
        case = copy.deepcopy(oracle.expected_cases()[0])
        baseline = oracle.configured_volume_mean_j_over_c(case)
        case["artificial_light_speed"] = 1250.0
        self.assertEqual(oracle.configured_volume_mean_j_over_c(case), baseline)
        case["species_mass"] = 7.0
        self.assertEqual(oracle.configured_volume_mean_j_over_c(case), baseline)
        with self.assertRaisesRegex(oracle.ContractError, "q/\\(mc\\) drifted"):
            oracle.validate_rendered_deck(case, oracle.render_oracle_deck(case))
        case["species_charge"] *= 2.0
        self.assertEqual(oracle.configured_volume_mean_j_over_c(case), 2.0 * baseline)

    def test_uniform_oracle_supports_arbitrary_positive_mass_with_consistent_q_over_mc(
        self,
    ) -> None:
        case = copy.deepcopy(oracle.expected_cases()[0])
        case["species_mass"] = 7.0
        case["species_charge"] = 7.0 * oracle.SPECIES_Q_OVER_MC
        case["species_charge_over_mass"] = oracle.SPECIES_Q_OVER_MC
        case["deposit_qscale"] = oracle.required_deposit_qscale(
            root_cell_volume=float(case["root_cell_volume"]),
            ppc=int(case["ppc"]),
            species_charge=float(case["species_charge"]),
        )
        validation = oracle.validate_rendered_deck(case, oracle.render_oracle_deck(case))
        self.assertEqual(validation["species_mass"], 7.0)
        self.assertAlmostEqual(
            validation["configured_volume_mean_j_over_c"], oracle.EXPECTED_J_OVER_C
        )

    def test_deck_rejects_species_velocity_overrides_that_change_current(self) -> None:
        case = _case("q043-current-oracle-d3-coarse-ppc1-single-cvr100")
        rendered = oracle.render_oracle_deck(case)
        nominal = oracle.parse_athinput_text(rendered)["particles"]
        failures = (
            ("zero", {"vx0": "0.0"}),
            ("different", {"vy0": str(float(nominal["cr_vy0"]) + 0.125)}),
        )
        for label, overrides in failures:
            with self.subTest(label=label):
                deck = _with_species0_overrides(rendered, overrides)
                with self.assertRaisesRegex(
                    oracle.ContractError, "species0 v[xy]0 must match nominal CR velocity"
                ):
                    oracle.validate_rendered_deck(case, deck)

    def test_deck_accepts_matching_species_velocity_overrides(self) -> None:
        case = _case("q043-current-oracle-d3-coarse-ppc1-single-cvr100")
        rendered = oracle.render_oracle_deck(case)
        nominal = oracle.parse_athinput_text(rendered)["particles"]
        deck = _with_species0_overrides(
            rendered,
            {
                "vx0": nominal["cr_vx0"],
                "vy0": nominal["cr_vy0"],
                "vz0": nominal["cr_vz0"],
            },
        )
        validation = oracle.validate_rendered_deck(case, deck)
        self.assertAlmostEqual(
            validation["configured_volume_mean_j_over_c"], oracle.EXPECTED_J_OVER_C
        )

    def test_deck_rejects_population_and_nominal_stream_drift(self) -> None:
        case = _case("q043-current-oracle-d1-coarse-ppc1-single-cvr100")
        rendered = oracle.render_oracle_deck(case)
        failures = (
            (
                "extra-species",
                rendered.replace("nspecies = 1", "nspecies = 2", 1),
                "requires exactly one species",
            ),
            (
                "reversed-nominal-stream",
                rendered.replace("cr_vx0 = 2.5", "cr_vx0 = -2.5", 1),
                "nominal CR velocity drifted",
            ),
        )
        for label, deck, message in failures:
            with self.subTest(label=label):
                with self.assertRaisesRegex(oracle.ContractError, message):
                    oracle.validate_rendered_deck(case, deck)

    def test_deck_and_helper_reject_fractional_ppc_without_truncation(self) -> None:
        case = oracle.expected_cases()[0]
        deck = oracle.render_oracle_deck(case).replace("ppc = 1.0", "ppc = 1.5", 1)
        with self.assertRaisesRegex(oracle.ContractError, "PPC must be a positive integer"):
            oracle.validate_rendered_deck(case, deck)
        with self.assertRaisesRegex(oracle.ContractError, "PPC must be a positive integer"):
            oracle.required_deposit_qscale(root_cell_volume=1.0, ppc=1.5)

    def test_cpp_source_binds_effective_species_velocity_to_nominal_stream(self) -> None:
        source = PGEN_SOURCE.read_text(encoding="utf-8")
        for axis in ("x", "y", "z"):
            with self.subTest(axis=axis):
                self.assertIn(
                    f'pin->GetReal("species0", "v{axis}0")',
                    source,
                )
                self.assertIn(
                    f'Q043VolumeAwareRequireClose("species0 v{axis}0", '
                    f"species_v{axis}, cr_v{axis});",
                    source,
                )
        self.assertIn(
            "species_charge, species_v_cr, root_cell_volume, b_g, k0",
            source,
        )

    def test_rendered_oracle_deck_satisfies_cpp_source_contract_end_to_end(self) -> None:
        case = oracle.expected_cases()[0]
        blocks = oracle.parse_athinput_text(oracle.render_oracle_deck(case))
        source = PGEN_SOURCE.read_text(encoding="utf-8")
        required_string_parameters = set(
            re.findall(
                r'Q043VolumeAwareRequireString\(\s*pin,\s*block,\s*"([^"]+)"',
                source,
                flags=re.DOTALL,
            )
        )
        required_string_parameters.update(
            re.findall(r'pin->GetString\(block,\s*"([^"]+)"', source)
        )
        self.assertFalse(required_string_parameters - set(blocks[oracle.PGEN_NAME]))
        self.assertEqual(
            blocks[oracle.PGEN_NAME]["q_over_mc_representation"],
            "species_charge_over_species_mass_equals_omega_over_b_g",
        )
        for index in range(1, len(oracle.FIELDS) + 1):
            self.assertEqual(
                blocks[f"output{index}"]["dcycle"], str(oracle.RAW_ORACLE_DCYCLE)
            )

        result = subprocess.run(
            [
                str(self.source_contract_binary),
                blocks[oracle.PGEN_NAME]["source_mode"],
                blocks[oracle.PGEN_NAME]["amplitude"],
                blocks["species0"]["mass"],
                blocks["species0"]["charge"],
                str(
                    float(blocks[oracle.PGEN_NAME]["omega"])
                    / float(blocks[oracle.PGEN_NAME]["b_g"])
                ),
            ],
            cwd=REPO_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            check=True,
        )
        self.assertEqual(result.stdout.strip(), "valid")

    def test_checked_in_decks_and_manifest_replay_exactly(self) -> None:
        manifest = oracle.validate_checked_in_decks()
        self.assertEqual(manifest["case_count"], oracle.EXPECTED_CASE_COUNT)
        self.assertFalse(manifest["launch_authorized"])
        self.assertFalse(manifest["scientific_claim_authorized"])
        self.assertFalse(manifest["publication_authorized"])
        self.assertFalse(manifest["artificial_c_in_formula"])
        self.assertEqual(manifest["required_output_cycle"], 1)
        self.assertEqual(manifest["output_dcycle"], 2)
        self.assertEqual(
            manifest["output_timing"],
            "cycle_zero_initialization_and_cycle_one_finalize_only",
        )
        self.assertEqual(
            manifest["initial_state"], "uniform_zero_perturbation_parallel_stream"
        )
        for record in manifest["cases"]:
            path = oracle.CHECKED_IN_DECK_ROOT / record["deck_path"]
            self.assertEqual(_sha256(path), record["deck_sha256"])
            self.assertAlmostEqual(
                record["configured_volume_mean_j_over_c"],
                oracle.EXPECTED_J_OVER_C,
            )

    def test_raw_single_and_multiaxis_cases_bind_provenance_and_pass_oracle(self) -> None:
        selected = (
            _case("q043-current-oracle-d1-coarse-ppc1-single-cvr100"),
            _case("q043-current-oracle-d2-coarse-ppc1-split_x2-cvr1000"),
            _case("q043-current-oracle-d3-fine-ppc4-split_xyz-cvr10000"),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for case in selected:
                with self.subTest(case=case["case_id"]):
                    paths = _write_raw_case(root, case)
                    report = oracle.analyze_raw_case(str(case["case_id"]), paths)
                    self.assertTrue(report["source_local_oracle_check_pass"])
                    self.assertFalse(report["passed"])
                    self.assertFalse(report["launch_authorized"])
                    self.assertFalse(report["scientific_claim_authorized"])
                    self.assertFalse(report["publication_authorized"])
                    self.assertEqual(
                        len(report["raw_provenance"]),
                        len(oracle.FIELDS) * int(case["mpi_ranks"]),
                    )
                    for binding in report["raw_provenance"]:
                        self.assertEqual(_sha256(Path(binding["path"])), binding["sha256"])
                    self.assertAlmostEqual(
                        report["measured_guide_projected_volume_mean_j_over_c"],
                        oracle.EXPECTED_J_OVER_C,
                        places=5,
                    )
                    self.assertLessEqual(
                        report["maximum_transverse_current"],
                        report["representation_derived_absolute_tolerance"],
                    )

    def test_complete_raw_matrix_verifies_all_axes_without_granting_authority(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            raw = _write_raw_matrix(Path(directory))
            report = oracle.analyze_raw_matrix(raw)
        self.assertTrue(report["source_local_oracle_check_pass"])
        self.assertFalse(report["passed"])
        self.assertEqual(report["case_count"], 132)
        self.assertEqual(report["dimensions_verified"], [1, 2, 3])
        self.assertEqual(report["resolutions_verified"], ["coarse", "fine"])
        self.assertEqual(report["ppc_verified"], [1, 4])
        self.assertEqual(report["decompositions_verified"], list(oracle.DECOMPOSITIONS))
        self.assertEqual(report["artificial_c_over_v_cr_verified"], [100, 1000, 10000])
        self.assertEqual(report["required_output_cycle"], 1)
        self.assertEqual(report["output_dcycle"], 2)
        self.assertFalse(report["artificial_c_in_formula"])
        self.assertFalse(report["launch_authorized"])
        self.assertFalse(report["scientific_claim_authorized"])
        self.assertFalse(report["publication_authorized"])

    def test_oracle_decks_fail_closed_if_output_dcycle_would_duplicate_cycle_one(self) -> None:
        case = oracle.expected_cases()[0]
        deck = oracle.render_oracle_deck(case).replace("dcycle = 2", "dcycle = 1")
        with self.assertRaisesRegex(oracle.ContractError, "output contract drifted"):
            oracle.validate_rendered_deck(case, deck)

    def test_real_runtime_output_bookkeeping_is_validated_and_normalized(
        self,
    ) -> None:
        case = _case("q043-current-oracle-d1-coarse-ppc1-single-cvr1000")
        self.assertEqual(
            oracle.expected_runtime_parameters(case)["particles"][
                "pic_boundary_conservation_ledger"
            ],
            "0",
        )
        self.assertEqual(
            oracle.expected_runtime_parameters(case)["particles"][
                "pic_allow_restart_injection_without_particle_section"
            ],
            "0",
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = _write_raw_case(root / "bookkeeping", case)
            report = oracle.analyze_raw_case(str(case["case_id"]), paths)
            self.assertTrue(report["source_local_oracle_check_pass"])
            self.assertEqual(
                report["cross_field_runtime_metadata"],
                oracle.RUNTIME_METADATA_CONTRACT,
            )

            drift: dict[str, dict[tuple[str, str], str]] = {"prtcl_jz": {}}
            drift["prtcl_jz"][("output4", "ghost_zones")] = "true"
            paths = _write_raw_case(root / "drift", case, parameter_overrides=drift)
            with self.assertRaisesRegex(oracle.ContractError, "immutable contract drifted"):
                oracle.analyze_raw_case(str(case["case_id"]), paths)

            enabled_ledger = {
                field: {("particles", "pic_boundary_conservation_ledger"): "1"}
                for field in oracle.FIELDS
            }
            paths = _write_raw_case(
                root / "enabled-ledger",
                case,
                parameter_overrides=enabled_ledger,
            )
            with self.assertRaisesRegex(
                oracle.ContractError,
                "authoritative deck and frozen default contract",
            ):
                oracle.analyze_raw_case(str(case["case_id"]), paths)

    def test_runtime_output_bookkeeping_rejects_malformed_or_impossible_states(
        self,
    ) -> None:
        case = _case("q043-current-oracle-d1-coarse-ppc1-single-cvr1000")
        failures = (
            (
                "malformed",
                {"prtcl_jx": {("output1", "file_number"): "not-an-integer"}},
                "file_number must be an integer",
            ),
            (
                "noncanonical",
                {"prtcl_jx": {("output1", "file_number"): "01"}},
                "canonical non-negative integer",
            ),
            (
                "impossible-progression",
                {"prtcl_jx": {("output3", "file_number"): "1"}},
                "sequential publication contract",
            ),
            (
                "invalid-last-time",
                {"prtcl_jx": {("output1", "last_time"): "0.0025"}},
                "cycle-cadence publication contract",
            ),
            (
                "common-mode-physics-drift",
                {
                    field: {("mhd", "gamma"): "9.0"}
                    for field in oracle.FIELDS
                },
                "authoritative deck and frozen default contract",
            ),
            (
                "common-mode-authority-drift",
                {
                    field: {
                        ("q043_bell_current_volume_aware", "launch_authorized"): "true"
                    }
                    for field in oracle.FIELDS
                },
                "authoritative deck and frozen default contract",
            ),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for label, overrides, message in failures:
                with self.subTest(label=label):
                    paths = _write_raw_case(
                        root / label, case, parameter_overrides=overrides
                    )
                    with self.assertRaisesRegex(oracle.ContractError, message):
                        oracle.analyze_raw_case(str(case["case_id"]), paths)

        parameters = oracle.parse_athinput_text(oracle.render_oracle_deck(case))
        parameters["output5"] = dict(parameters["output4"])
        with self.assertRaisesRegex(oracle.ContractError, "output block inventory"):
            oracle._normalized_runtime_parameters(
                parameters, case=case, field="prtcl_rho"
            )

    def test_runtime_accepts_exact_explicit_multidimensional_species_velocity(
        self,
    ) -> None:
        case = _case("q043-current-oracle-d2-coarse-ppc1-single-cvr1000")
        particles = oracle.parse_athinput_text(oracle.render_oracle_deck(case))[
            "particles"
        ]
        exact_velocity = {
            field: {
                ("species0", f"v{axis}0"): particles[f"cr_v{axis}0"]
                for axis in ("x", "y", "z")
            }
            for field in oracle.FIELDS
        }
        with tempfile.TemporaryDirectory() as directory:
            paths = _write_raw_case(
                Path(directory),
                case,
                parameter_overrides=exact_velocity,
            )
            report = oracle.analyze_raw_case(str(case["case_id"]), paths)
            self.assertTrue(report["source_local_oracle_check_pass"])

    def test_transverse_current_and_spatial_nonuniformity_fail_closed(self) -> None:
        case = _case("q043-current-oracle-d1-coarse-ppc1-single-cvr1000")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            transverse = _write_raw_case(
                root / "transverse",
                case,
                mutations={("prtcl_jy", 0): (0, 0.1)},
            )
            with self.assertRaisesRegex(oracle.ContractError, "transverse deposited"):
                oracle.analyze_raw_case(str(case["case_id"]), transverse)

            nonuniform = _write_raw_case(
                root / "nonuniform",
                case,
                mutations={("prtcl_jx", 0): (0, 0.001)},
            )
            with self.assertRaisesRegex(oracle.ContractError, "spatial nonuniformity"):
                oracle.analyze_raw_case(str(case["case_id"]), nonuniform)

            rho = _write_raw_case(
                root / "rho",
                case,
                mutations={("prtcl_rho", 0): (0, 0.001)},
            )
            with self.assertRaisesRegex(oracle.ContractError, "rho spatial nonuniformity"):
                oracle.analyze_raw_case(str(case["case_id"]), rho)

    def test_matrix_completeness_shards_and_raw_reuse_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = _write_raw_matrix(root)
            raw.pop(next(iter(raw)))
            with self.assertRaisesRegex(oracle.ContractError, "matrix is incomplete"):
                oracle.analyze_raw_matrix(raw)

            split = _case("q043-current-oracle-d3-coarse-ppc1-split_xyz-cvr1000")
            split_paths = _write_raw_case(root / "split", split)
            split_paths["prtcl_jx"] = split_paths["prtcl_jx"][:1]
            with self.assertRaisesRegex(oracle.ContractError, "shard count"):
                oracle.analyze_raw_case(str(split["case_id"]), split_paths)

            raw = _write_raw_matrix(root / "reuse")
            ids = list(raw)
            raw[ids[1]] = raw[ids[0]]
            with self.assertRaises(oracle.ContractError):
                oracle.analyze_raw_matrix(raw)

    def test_readiness_record_binds_additive_non_authorizing_oracle(self) -> None:
        record = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertEqual(record["campaign_id"], oracle.CAMPAIGN_ID)
        self.assertEqual(record["qualification_effect"], oracle.QUALIFICATION_EFFECT)
        self.assertFalse(record["authority"]["launch_authorized"])
        self.assertFalse(record["authority"]["scientific_claim_authorized"])
        self.assertFalse(record["authority"]["publication_authorized"])
        self.assertEqual(record["matrix_contract"]["case_count"], 132)
        self.assertEqual(
            record["matrix_contract"]["decompositions"],
            list(oracle.DECOMPOSITIONS),
        )
        self.assertEqual(
            record["matrix_contract"]["decompositions_by_dimension"],
            {
                str(dimension): list(oracle.DECOMPOSITIONS_BY_DIMENSION[dimension])
                for dimension in oracle.DIMENSIONS
            },
        )
        self.assertEqual(
            record["normalization_contract"]["formula"],
            "PPC*deposit_qscale*species_charge*v_CR/V_root_cell=2*B_g*k0",
        )
        self.assertEqual(
            record["normalization_contract"]["species_population_contract"],
            "exactly_one_species",
        )
        self.assertIn(
            "effective_species0_vx0_vy0_vz0_must_match_nominal",
            record["normalization_contract"]["species_velocity_contract"],
        )
        self.assertIn(
            "positive_integer_only",
            record["normalization_contract"]["ppc_contract"],
        )
        self.assertEqual(record["raw_output_oracle"]["cycle"], 1)
        self.assertEqual(record["raw_output_oracle"]["output_dcycle"], 2)
        self.assertEqual(
            record["raw_output_oracle"]["cross_field_runtime_metadata"],
            oracle.RUNTIME_METADATA_CONTRACT,
        )
        self.assertEqual(
            record["raw_output_oracle"]["initial_state"],
            "uniform_zero_perturbation_parallel_stream",
        )
        self.assertFalse(
            record["normalization_contract"]["artificial_light_speed_in_formula"]
        )
        self.assertEqual(
            _sha256(oracle.CHECKED_IN_MANIFEST),
            record["artifact_bindings"][
                "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle/deck_manifest.json"
            ],
        )
        for relative, digest in record["artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), digest)


if __name__ == "__main__":
    unittest.main()
