#!/usr/bin/env python3
"""Focused tests for exact source-local Q-023 Bell deck materialization."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
import uuid

from tst.publication import materialize_q023_paper_bell_linear_variants as materializer


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _parse_deck_text(text: str) -> dict[str, dict[str, str]]:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "variant.athinput"
        path.write_text(text, encoding="utf-8")
        return materializer.bell.parse_athinput(path)


class Q023PaperBellLinearMaterializerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.request = materializer.build_materialization_request()
        (
            cls.manifest,
            cls.decks,
            cls.request_text,
        ) = materializer.build_materialization_manifest(cls.request)

    def _selector(
        self,
        dimension: int,
        epsilon: float,
        resolution_scale: float,
        timestep_scale: float,
        ppc: int,
    ) -> dict[str, object]:
        return next(
            item
            for item in self.request["variants"]
            if (
                item["dimension"],
                item["epsilon"],
                item["resolution_scale"],
                item["timestep_scale"],
                item["ppc"],
            )
            == (dimension, epsilon, resolution_scale, timestep_scale, ppc)
        )

    def test_complete_fixed_grid_manifest_and_decks_are_deterministic(self) -> None:
        selectors = materializer.validate_materialization_request(self.request)
        self.assertEqual(
            len(selectors), materializer.EXPECTED_VARIANT_COUNT
        )
        self.assertEqual(len(self.manifest["variants"]), 405)
        self.assertEqual(len(self.decks), 405)
        self.assertFalse(self.manifest["qualifying_evidence"])
        self.assertEqual(
            self.manifest["qualification_effect"],
            materializer.QUALIFICATION_EFFECT,
        )
        self.assertEqual(
            self.manifest["section52_qualification"], "not_claimed"
        )
        self.assertEqual(
            self.manifest["loading_policy"], materializer.LOADING_POLICY
        )
        self.assertIn("freeze_still_open", self.manifest["timestep_policy"])
        self.assertEqual(
            self.manifest["materialization_request_sha256"],
            _sha256_bytes(self.request_text.encode("utf-8")),
        )

        manifest, decks, request_text = (
            materializer.build_materialization_manifest(self.request)
        )
        self.assertEqual(manifest, self.manifest)
        self.assertEqual(decks, self.decks)
        self.assertEqual(request_text, self.request_text)

        for variant in self.manifest["variants"]:
            self.assertEqual(variant["cr_distribution"], "center")
            self.assertFalse(variant["qualifying_evidence"])
            self.assertEqual(variant["section52_qualification"], "not_claimed")
            self.assertAlmostEqual(
                variant["ppc"] * variant["deposit_qscale"],
                2.0e9,
            )
            speed = math.sqrt(
                sum(value * value for value in variant["stream_velocity"])
            )
            self.assertAlmostEqual(speed, 1.0 / variant["epsilon"])
            self.assertAlmostEqual(
                variant["pic_cr_light_speed"], 1000.0 * speed
            )
            self.assertEqual(
                _sha256_bytes(self.decks[variant["deck_path"]].encode("utf-8")),
                variant["deck_sha256"],
            )

    def test_rendered_variants_scale_epsilon_resolution_timestep_and_ppc(
        self,
    ) -> None:
        one_d, one_d_text = materializer.render_variant_deck(
            self._selector(1, 0.1, 2.0, 0.5, 32)
        )
        blocks = _parse_deck_text(one_d_text)
        self.assertEqual(blocks["mesh"]["nx1"], "64")
        self.assertEqual(blocks["mesh"]["nx2"], "4")
        self.assertEqual(blocks["meshblock"]["nx1"], "64")
        self.assertEqual(blocks["meshblock"]["nx2"], "4")
        self.assertEqual(blocks["time"]["cfl_number"], "0.05")
        self.assertEqual(blocks["particles"]["ppc"], "32.0")
        self.assertEqual(blocks["particles"]["deposit_qscale"], "62500000.0")
        self.assertEqual(blocks["particles"]["cr_vx0"], "10.0")
        self.assertEqual(blocks["particles"]["pic_cr_light_speed"], "10000.0")
        self.assertEqual(blocks["particles"]["cr_distribution"], "center")
        self.assertEqual(
            blocks["q023_paper_bell_linear"]["timestep"],
            "open_clean_candidate_timestep_freeze",
        )
        self.assertEqual(one_d["deposit_qscale"], 62500000.0)

        _, three_d_text = materializer.render_variant_deck(
            self._selector(3, 0.8, 0.5, 2.0, 8)
        )
        blocks = _parse_deck_text(three_d_text)
        self.assertEqual(
            [blocks["mesh"][f"nx{axis}"] for axis in (1, 2, 3)],
            ["64", "32", "16"],
        )
        self.assertEqual(
            [blocks["meshblock"][f"nx{axis}"] for axis in (1, 2, 3)],
            ["16", "16", "16"],
        )
        self.assertEqual(blocks["time"]["cfl_number"], "0.2")
        self.assertEqual(blocks["particles"]["deposit_qscale"], "250000000.0")
        self.assertEqual(blocks["particles"]["pic_cr_light_speed"], "1250.0")
        self.assertEqual(blocks["q023_paper_bell_linear"]["epsilon"], "0.8")

    def test_materialization_writes_only_below_safe_new_codex_root(self) -> None:
        output_root = (
            materializer.AUTHORIZED_OUTPUT_PARENT
            / f"q023-test-materialized-{uuid.uuid4().hex}"
        )
        try:
            manifest = materializer.materialize_variant_decks(
                output_root, self.request
            )
            self.assertEqual(manifest, self.manifest)
            self.assertEqual(
                (output_root / "materialization_request.json").read_text(
                    encoding="utf-8"
                ),
                self.request_text,
            )
            written_manifest = json.loads(
                (output_root / "materialization_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(written_manifest, self.manifest)
            files = [path for path in output_root.rglob("*") if path.is_file()]
            self.assertEqual(len(files), materializer.EXPECTED_VARIANT_COUNT + 2)
            selected = self.manifest["variants"][0]
            self.assertEqual(
                _sha256_bytes((output_root / selected["deck_path"]).read_bytes()),
                selected["deck_sha256"],
            )
            with self.assertRaisesRegex(
                materializer.ContractError, "already exists"
            ):
                materializer.materialize_variant_decks(output_root, self.request)
        finally:
            shutil.rmtree(output_root, ignore_errors=True)

    def test_unknown_duplicate_incomplete_and_qualification_requests_fail_closed(
        self,
    ) -> None:
        request = copy.deepcopy(self.request)
        request["variants"].append(copy.deepcopy(request["variants"][0]))
        with self.assertRaisesRegex(materializer.ContractError, "duplicate"):
            materializer.validate_materialization_request(request)

        request = copy.deepcopy(self.request)
        request["variants"].pop()
        with self.assertRaisesRegex(materializer.ContractError, "incomplete"):
            materializer.validate_materialization_request(request)

        request = copy.deepcopy(self.request)
        request["runtime_artifact_root"] = "/tmp/not-admitted"
        with self.assertRaisesRegex(materializer.ContractError, "fields"):
            materializer.validate_materialization_request(request)

        request = copy.deepcopy(self.request)
        request["variants"][0]["runtime_artifact_path"] = "/tmp/not-admitted"
        with self.assertRaisesRegex(materializer.ContractError, "fields"):
            materializer.validate_materialization_request(request)

        request = copy.deepcopy(self.request)
        request["section52_qualification"] = "claimed"
        with self.assertRaisesRegex(materializer.ContractError, "refuses"):
            materializer.validate_materialization_request(request)

        request = copy.deepcopy(self.request)
        request["qualification_manifest"] = {"status": "claimed"}
        with self.assertRaisesRegex(materializer.ContractError, "refuses"):
            materializer.validate_materialization_request(request)

        request = copy.deepcopy(self.request)
        request["variants"][0]["qualification"] = "claimed"
        with self.assertRaisesRegex(materializer.ContractError, "refuses"):
            materializer.validate_materialization_request(request)

        request = copy.deepcopy(self.request)
        request["loading_policy"]["cr_distribution"] = "random"
        with self.assertRaisesRegex(materializer.ContractError, "centered-loading"):
            materializer.validate_materialization_request(request)

    def test_unsafe_output_roots_fail_closed(self) -> None:
        with self.assertRaisesRegex(materializer.ContractError, "absolute"):
            materializer.materialize_variant_decks(
                Path("relative-materialized-root"), self.request
            )

        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(materializer.ContractError, "outside"):
                materializer.materialize_variant_decks(
                    Path(directory) / "materialized", self.request
                )

        with tempfile.TemporaryDirectory(
            dir=materializer.AUTHORIZED_OUTPUT_PARENT
        ) as directory:
            with self.assertRaisesRegex(materializer.ContractError, "direct"):
                materializer.materialize_variant_decks(
                    Path(directory) / "materialized", self.request
                )

        with tempfile.TemporaryDirectory(
            dir=materializer.AUTHORIZED_OUTPUT_PARENT
        ) as directory:
            directory_path = Path(directory)
            actual = directory_path / "actual"
            actual.mkdir()
            alias = directory_path / "alias"
            alias.symlink_to(actual, target_is_directory=True)
            with self.assertRaisesRegex(materializer.ContractError, "direct"):
                materializer.materialize_variant_decks(
                    alias / "materialized", self.request
                )

    def test_reviewed_source_deck_digest_drift_fails_closed(self) -> None:
        with patch.dict(
            materializer.SOURCE_DECK_SHA256, {2: "0" * 64}
        ):
            with self.assertRaisesRegex(
                materializer.ContractError, "reviewed 2D source deck digest"
            ):
                materializer.build_materialization_manifest(self.request)


if __name__ == "__main__":
    unittest.main()
