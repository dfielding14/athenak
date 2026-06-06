#!/usr/bin/env python3
"""Focused adversarial tests for the additive Q022 central-reference maps."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "tst/publication/q022_central_reference_mapping_successor_v1.py"
SPEC = importlib.util.spec_from_file_location("q022_central_reference_mapping_successor_v1", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
q022 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(q022)


class Q022CentralReferenceMappingSuccessorTests(unittest.TestCase):
    @staticmethod
    def _load(path: Path) -> dict[str, object]:
        return json.loads(path.read_text(encoding="utf-8"))

    def test_complete_successor_validates(self) -> None:
        q022.validate_all()

    def test_artificial_C_cannot_return_to_Bell_target(self) -> None:
        record = self._load(q022.BELL_MAP_PATH)
        record["normalization_map"]["artificial_C_multiplier_allowed"] = True
        with self.assertRaisesRegex(q022.ContractError, "legacy Bell normalization"):
            q022.validate_bell_map(record)

    def test_Bell_map_cannot_claim_standalone_nonlinear_reference_dataset(self) -> None:
        record = self._load(q022.BELL_MAP_PATH)
        record["parameter_overlap"]["standalone_nonlinear_reference_dataset_available"] = True
        with self.assertRaisesRegex(q022.ContractError, "standalone nonlinear dataset"):
            q022.validate_bell_map(record)

    def test_Sun_Bai_pressure_cannot_be_relabelled_explicit(self) -> None:
        record = self._load(q022.SHOCK_MAP_PATH)
        record["parameter_overlap"]["not_explicit_in_Sun_Bai_Section_5_4"] = {}
        with self.assertRaisesRegex(q022.ContractError, "pressure uncertainty"):
            q022.validate_shock_map(record)

    def test_numeric_tolerance_cannot_be_invented(self) -> None:
        record = self._load(q022.SHOCK_MAP_PATH)
        record["comparison_tolerances"] = [{"observable": "Bmax/B0", "value": 0.1}]
        with self.assertRaisesRegex(q022.ContractError, "numeric tolerances"):
            q022.validate_shock_map(record)

    def test_page_provenance_must_address_bound_PDF(self) -> None:
        record = self._load(q022.BELL_MAP_PATH)
        record["provenance_entries"][0]["physical_pdf_page"] = 0
        with self.assertRaisesRegex(q022.ContractError, "physical PDF page"):
            q022.validate_bell_map(record)

    def test_reference_PDF_hash_must_match(self) -> None:
        record = self._load(q022.SHOCK_MAP_PATH)
        record["reference_artifacts"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(q022.ContractError, "reference PDF drifted"):
            q022.validate_shock_map(record)

    def test_authority_must_remain_false(self) -> None:
        record = self._load(q022.LEDGER_PATH)
        record["authority"]["execution_authorized"] = True
        with self.assertRaisesRegex(q022.ContractError, "authority must remain false"):
            q022.validate_ledger(record)

    def test_ledger_inventory_is_fail_closed(self) -> None:
        record = self._load(q022.LEDGER_PATH)
        record["entries"] = record["entries"][1:]
        with self.assertRaisesRegex(q022.ContractError, "ledger inventory drifted"):
            q022.validate_ledger(record)

    def test_readiness_binds_historical_and_current_contracts(self) -> None:
        record = self._load(q022.READINESS_PATH)
        q022.validate_readiness(record)
        mutated = copy.deepcopy(record)
        mutated["historical_q022_bindings"][0]["sha256"] = "f" * 64
        with self.assertRaisesRegex(q022.ContractError, "bound artifact drifted"):
            q022.validate_readiness(mutated)
        missing = copy.deepcopy(record)
        missing["historical_q022_bindings"] = missing["historical_q022_bindings"][1:]
        with self.assertRaisesRegex(q022.ContractError, "historical Q022 inventory"):
            q022.validate_readiness(missing)


if __name__ == "__main__":
    unittest.main()
