#!/usr/bin/env python3
"""Focused tests for the corrected Bell deposited-current raw-output oracle."""

from __future__ import annotations

import copy
import hashlib
import inspect
import json
import math
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

from tst.publication import q023_bell_joverc_deposited_current_oracle as oracle


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q023_bell_joverc_deposited_current_oracle_source_local_2026-06-06.json"
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


def _binary_payload(case: dict[str, object], field: str, blocks: list[bytes]) -> bytes:
    parameter_header = oracle.render_oracle_deck(case).encode("utf-8")
    return (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        b"  time=0.0\n"
        b"  cycle=0\n"
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
) -> dict[str, tuple[Path, ...]]:
    mutations = mutations or {}
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


class Q023BellJOverCDepositedCurrentOracleTests(unittest.TestCase):
    def test_exact_deck_matrix_derives_qscale_from_root_volume_and_ppc(self) -> None:
        cases = oracle.expected_cases()
        self.assertEqual(len(cases), 24)
        self.assertEqual({case["dimension"] for case in cases}, {1, 2, 3})
        self.assertEqual({case["resolution"] for case in cases}, {"coarse", "fine"})
        self.assertEqual({case["ppc"] for case in cases}, {1, 4})
        self.assertEqual({case["decomposition"] for case in cases}, {"single", "split"})
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
                / (int(case["ppc"]) * oracle.CHARGE_Q_OVER_MC * oracle.STREAM_SPEED)
            )
            self.assertAlmostEqual(float(case["deposit_qscale"]), expected)
            self.assertAlmostEqual(
                oracle.configured_volume_mean_j_over_c(case),
                oracle.EXPECTED_J_OVER_C,
            )
        for pair in grouped.values():
            self.assertEqual(len(pair), 2)
            self.assertEqual(pair[0]["deposit_qscale"], pair[1]["deposit_qscale"])

    def test_artificial_c_is_absent_from_configured_current_formula(self) -> None:
        signature = inspect.signature(oracle.required_deposit_qscale)
        self.assertNotIn("artificial_light_speed", signature.parameters)
        case = copy.deepcopy(oracle.expected_cases()[0])
        baseline = oracle.configured_volume_mean_j_over_c(case)
        case["artificial_light_speed"] = 1250.0
        self.assertEqual(oracle.configured_volume_mean_j_over_c(case), baseline)

    def test_checked_in_decks_and_manifest_replay_exactly(self) -> None:
        manifest = oracle.validate_checked_in_decks()
        self.assertEqual(manifest["case_count"], oracle.EXPECTED_CASE_COUNT)
        self.assertFalse(manifest["launch_authorized"])
        self.assertFalse(manifest["scientific_claim_authorized"])
        self.assertFalse(manifest["publication_authorized"])
        self.assertFalse(manifest["artificial_c_in_formula"])
        for record in manifest["cases"]:
            path = oracle.CHECKED_IN_DECK_ROOT / record["deck_path"]
            self.assertEqual(_sha256(path), record["deck_sha256"])
            self.assertAlmostEqual(
                record["configured_volume_mean_j_over_c"],
                oracle.EXPECTED_J_OVER_C,
            )

    def test_raw_single_and_split_cases_bind_provenance_and_pass_oracle(self) -> None:
        selected = (
            _case("q023-joverc-current-oracle-d1-coarse-ppc1-single"),
            _case("q023-joverc-current-oracle-d3-fine-ppc4-split"),
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
        self.assertEqual(report["case_count"], 24)
        self.assertEqual(report["dimensions_verified"], [1, 2, 3])
        self.assertEqual(report["resolutions_verified"], ["coarse", "fine"])
        self.assertEqual(report["ppc_verified"], [1, 4])
        self.assertEqual(report["decompositions_verified"], ["single", "split"])
        self.assertFalse(report["artificial_c_in_formula"])
        self.assertFalse(report["launch_authorized"])
        self.assertFalse(report["scientific_claim_authorized"])
        self.assertFalse(report["publication_authorized"])

    def test_transverse_current_and_spatial_nonuniformity_fail_closed(self) -> None:
        case = _case("q023-joverc-current-oracle-d1-coarse-ppc1-single")
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

            split = _case("q023-joverc-current-oracle-d2-coarse-ppc1-split")
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
        self.assertEqual(record["matrix_contract"]["case_count"], 24)
        self.assertEqual(
            record["normalization_contract"]["formula"],
            "PPC*deposit_qscale*species_charge*v_CR/V_root_cell=2*B_g*k0",
        )
        self.assertFalse(
            record["normalization_contract"]["artificial_light_speed_in_formula"]
        )
        self.assertEqual(
            _sha256(oracle.CHECKED_IN_MANIFEST),
            record["artifact_bindings"][
                "inputs/tests/q023_bell_joverc_deposited_current_oracle/deck_manifest.json"
            ],
        )
        for relative, digest in record["artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), digest)


if __name__ == "__main__":
    unittest.main()
