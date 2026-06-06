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


def _block_payload(
    case: dict[str, object],
    field: str,
    logical_x1: int,
    *,
    values: np.ndarray,
) -> bytes:
    nx1, nx2, nx3 = (int(value) for value in case["meshblock_nx"])
    bounds = tuple(tuple(float(value) for value in item) for item in case["bounds"])
    root_blocks_x1 = int(case["global_nx"][0]) // nx1
    x1_width = (bounds[0][1] - bounds[0][0]) / root_blocks_x1
    geometry = (
        bounds[0][0] + logical_x1 * x1_width,
        bounds[0][0] + (logical_x1 + 1) * x1_width,
        bounds[1][0],
        bounds[1][1],
        bounds[2][0],
        bounds[2][1],
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
        logical_x1,
        0,
        0,
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
    pending = dict(overrides)
    rendered = []
    block: str | None = None
    for line in oracle.render_oracle_deck(case).splitlines():
        if line.startswith("<") and line.endswith(">"):
            if block is not None:
                for (candidate, key), value in list(pending.items()):
                    if candidate == block:
                        rendered.append(f"{key} = {value}")
                        del pending[(candidate, key)]
            block = line[1:-1]
            rendered.append(line)
            continue
        if block is not None and "=" in line:
            key = line.split("=", 1)[0].strip()
            override = pending.pop((block, key), None)
            if override is not None:
                rendered.append(f"{key} = {override}")
                continue
        rendered.append(line)
    if block is not None:
        for (candidate, key), value in list(pending.items()):
            if candidate == block:
                rendered.append(f"{key} = {value}")
                del pending[(candidate, key)]
    if pending:
        raise AssertionError(f"unknown runtime parameter override(s): {sorted(pending)}")
    return "\n".join(rendered)


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
    for field in oracle.FIELDS:
        field_paths = []
        for shard in range(shard_count):
            values = np.full(shape, _field_value(case, field), dtype=np.float64)
            mutation = mutations.get((field, shard))
            if mutation is not None:
                flat_index, delta = mutation
                values.reshape(-1)[flat_index] += delta
            payload = _binary_payload(
                case,
                field,
                [_block_payload(case, field, shard, values=values)],
                parameter_overrides=parameter_overrides.get(field),
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
        self.assertEqual(len(cases), 72)
        self.assertEqual({case["dimension"] for case in cases}, {1, 2, 3})
        self.assertEqual({case["resolution"] for case in cases}, {"coarse", "fine"})
        self.assertEqual({case["ppc"] for case in cases}, {1, 4})
        self.assertEqual({case["decomposition"] for case in cases}, {"single", "split"})
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
        for pair in grouped.values():
            self.assertEqual(len(pair), 6)
            self.assertEqual({item["deposit_qscale"] for item in pair}, {pair[0]["deposit_qscale"]})

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

    def test_raw_single_and_split_cases_bind_provenance_and_pass_oracle(self) -> None:
        selected = (
            _case("q043-current-oracle-d1-coarse-ppc1-single-cvr100"),
            _case("q043-current-oracle-d3-fine-ppc4-split-cvr10000"),
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
        self.assertEqual(report["case_count"], 72)
        self.assertEqual(report["dimensions_verified"], [1, 2, 3])
        self.assertEqual(report["resolutions_verified"], ["coarse", "fine"])
        self.assertEqual(report["ppc_verified"], [1, 4])
        self.assertEqual(report["decompositions_verified"], ["single", "split"])
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

    def test_real_runtime_output_bookkeeping_is_normalized_but_other_drift_fails(
        self,
    ) -> None:
        case = _case("q043-current-oracle-d1-coarse-ppc1-single-cvr1000")
        bookkeeping = {}
        for field_index, field in enumerate(oracle.FIELDS):
            bookkeeping[field] = {
                (f"output{output_index}", "file_number"): str(
                    1 + int(output_index <= field_index + 1)
                )
                for output_index in range(1, len(oracle.FIELDS) + 1)
            }
            bookkeeping[field].update(
                {
                    (f"output{output_index}", "last_time"): str(field_index)
                    for output_index in range(1, len(oracle.FIELDS) + 1)
                }
            )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = _write_raw_case(root / "bookkeeping", case, parameter_overrides=bookkeeping)
            report = oracle.analyze_raw_case(str(case["case_id"]), paths)
            self.assertTrue(report["source_local_oracle_check_pass"])
            self.assertEqual(
                report["cross_field_runtime_metadata"],
                "strict_after_normalizing_only_numbered_output_file_number_and_last_time",
            )

            drift = copy.deepcopy(bookkeeping)
            drift["prtcl_jz"][("output4", "ghost_zones")] = "true"
            paths = _write_raw_case(root / "drift", case, parameter_overrides=drift)
            with self.assertRaisesRegex(oracle.ContractError, "raw field metadata disagrees"):
                oracle.analyze_raw_case(str(case["case_id"]), paths)

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

            split = _case("q043-current-oracle-d2-coarse-ppc1-split-cvr1000")
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
        self.assertEqual(record["matrix_contract"]["case_count"], 72)
        self.assertEqual(
            record["normalization_contract"]["formula"],
            "PPC*deposit_qscale*species_charge*v_CR/V_root_cell=2*B_g*k0",
        )
        self.assertEqual(record["raw_output_oracle"]["cycle"], 1)
        self.assertEqual(record["raw_output_oracle"]["output_dcycle"], 2)
        self.assertEqual(
            record["raw_output_oracle"]["cross_field_runtime_metadata"],
            "strict_after_normalizing_only_numbered_output_file_number_and_last_time",
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
