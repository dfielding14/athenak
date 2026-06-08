#!/usr/bin/env python3
"""Focused tests for the registered Q023 retirement driver."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile

import pytest

from tst.publication import q023_registered_retirement_driver_v1 as driver


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")
    path.chmod(0o444)


def test_retired_active_pair_requires_exact_mirrors_and_lifecycle() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        orion = root / "orion"
        project = root / "project"
        version = "a" * 64
        policy = {
            "registered_science_slices": [],
            "olcf_side_storage": {"installed_control_plane_version": version},
        }
        for target in (orion, project):
            _write_json(target / "policy/storage_policy.json", policy)
            _write_json(target / "policy/active_promotion.json", {"active": True})
        result = driver.validate_retired_active_pair(
            policy,
            control_plane_version=version,
            pic_root=orion,
            project_home_root=project,
        )
        assert all(driver.SHA256_PATTERN.fullmatch(value) for value in result)
        drifted = {**policy, "registered_science_slices": [{"unexpected": True}]}
        (project / "policy/storage_policy.json").chmod(0o644)
        _write_json(project / "policy/storage_policy.json", drifted)
        with pytest.raises(driver.RetirementError, match="mirrors differ"):
            driver.validate_retired_active_pair(
                policy,
                control_plane_version=version,
                pic_root=orion,
                project_home_root=project,
            )


def test_stable_json_rejects_writable_or_nonobject_files() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        writable = root / "writable.json"
        writable.write_text("{}\n", encoding="utf-8")
        with pytest.raises(driver.RetirementError, match="read-only"):
            driver._stable_read_only_json(writable)
        array = root / "array.json"
        array.write_text("[]\n", encoding="utf-8")
        array.chmod(0o444)
        with pytest.raises(driver.RetirementError, match="JSON object"):
            driver._stable_read_only_json(array)
