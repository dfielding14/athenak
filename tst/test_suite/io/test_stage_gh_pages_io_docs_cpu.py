"""Focused safety tests for the detached gh-pages IO documentation helper."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import stat
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
HELPER_PATH = REPO_ROOT / "scripts" / "stage_gh_pages_io_docs.py"
SPEC = importlib.util.spec_from_file_location("stage_gh_pages_io_docs", HELPER_PATH)
assert SPEC is not None and SPEC.loader is not None
helper = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = helper
SPEC.loader.exec_module(helper)


def _run(repo: Path, *args: str) -> str:
    return subprocess.run(
        [*args],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    ).stdout.strip()


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _blob(repo: Path, path: str) -> str:
    return _run(repo, "git", "rev-parse", f"HEAD:{path}")


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


@pytest.fixture
def pages_fixture(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    _run(source, "git", "init", "-b", "main")
    _run(source, "git", "config", "user.name", "RCP-07 Test")
    _run(source, "git", "config", "user.email", "rcp07@example.invalid")
    baseline = {
        ".gitignore": "ignored.rcp07\n",
        "docs/source/index.md": "# Protected Index\n",
        "docs/source/page.md": "# Page\n\n## Replace\nold\n\n## Next\nkeep\n",
        "docs/source/insert.md": "# Insert\n\nbefore\n\n## Anchor\nkeep\n",
        "docs/source/whole.md": "old whole\n",
        "docs/source/unrelated.md": "untouched\n",
    }
    for path, text in baseline.items():
        _write(source / path, text)
    _run(source, "git", "add", ".")
    _run(source, "git", "commit", "-m", "baseline")
    head = _run(source, "git", "rev-parse", "HEAD")
    _run(source, "git", "update-ref", "refs/remotes/origin/gh-pages", head)

    bundle = source / "bundle"
    payloads = {
        "fragments/page.md": "## Replace\nnew-feature\n",
        "fragments/insert.md": "## Inserted\ninserted\n",
        "overlay/whole.md": "new whole\n",
        "overlay/add.md": "new add\n",
    }
    for path, text in payloads.items():
        _write(bundle / path, text)
    manifest = {
        "schema_version": 1,
        "bundle": "test",
        "expected_baseline": head,
        "allowlist": [
            "docs/source/page.md",
            "docs/source/insert.md",
            "docs/source/whole.md",
            "docs/source/add.md",
        ],
        "protected_blobs": {
            "docs/source/index.md": _blob(source, "docs/source/index.md"),
        },
        "targets": [
            {
                "path": "docs/source/page.md",
                "baseline_blob": _blob(source, "docs/source/page.md"),
                "operations": [
                    {
                        "kind": "replace_section",
                        "marker": "page",
                        "source": "fragments/page.md",
                        "source_sha256": _sha256(payloads["fragments/page.md"].encode()),
                        "start_anchor": "## Replace\n",
                        "end_anchor": "\n## Next\n",
                        "baseline_section_sha256": _sha256(b"## Replace\nold\n"),
                    }
                ],
            },
            {
                "path": "docs/source/insert.md",
                "baseline_blob": _blob(source, "docs/source/insert.md"),
                "operations": [
                    {
                        "kind": "insert_before",
                        "marker": "insert",
                        "source": "fragments/insert.md",
                        "source_sha256": _sha256(
                            payloads["fragments/insert.md"].encode()
                        ),
                        "anchor": "\n## Anchor\n",
                        "anchor_sha256": _sha256(b"\n## Anchor\n"),
                    }
                ],
            },
            {
                "path": "docs/source/whole.md",
                "baseline_blob": _blob(source, "docs/source/whole.md"),
                "operations": [
                    {
                        "kind": "replace_file",
                        "source": "overlay/whole.md",
                        "source_sha256": _sha256(payloads["overlay/whole.md"].encode()),
                    }
                ],
            },
            {
                "path": "docs/source/add.md",
                "expected_absent": True,
                "operations": [
                    {
                        "kind": "add_file",
                        "source": "overlay/add.md",
                        "source_sha256": _sha256(payloads["overlay/add.md"].encode()),
                    }
                ],
            },
        ],
        "contradiction_checks": [
            {
                "pattern": "stale-guidance",
                "description": "stale test guidance",
                "maximum": 0,
            },
            {
                "pattern": "new-feature",
                "description": "replacement test guidance",
                "minimum": 1,
            },
        ],
    }
    manifest_path = bundle / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    target = tmp_path / "target"
    _run(source, "git", "worktree", "add", "--detach", str(target), head)
    return source, target, bundle, manifest_path, manifest


def _sources(manifest, bundle):
    helper._validate_manifest(manifest, bundle)
    return helper._load_sources(manifest, bundle)


def _mutate_once_during_packet(monkeypatch, path, text):
    original = helper._write_packet_file
    mutated = False

    def write_packet_file(*args, **kwargs):
        nonlocal mutated
        if not mutated:
            mutated = True
            _write(path, text)
        return original(*args, **kwargs)

    monkeypatch.setattr(helper, "_write_packet_file", write_packet_file)


def test_clean_stage_verify_and_preserve_unrelated_content(pages_fixture):
    source, target, _, manifest_path, manifest = pages_fixture
    helper.execute(target, source_repo=source, manifest_path=manifest_path)

    page = (target / "docs/source/page.md").read_text()
    insertion = (target / "docs/source/insert.md").read_text()
    assert f"<!-- BEGIN {helper.MARKER_PREFIX}: page -->" in page
    assert f"<!-- END {helper.MARKER_PREFIX}: page -->\n## Next" in page
    assert "## Next\nkeep" in page
    assert f"<!-- BEGIN {helper.MARKER_PREFIX}: insert -->" in insertion
    assert insertion.index("## Inserted") < insertion.index("## Anchor")
    assert (target / "docs/source/unrelated.md").read_text() == "untouched\n"
    assert (target / "docs/source/add.md").read_text() == "new add\n"
    assert helper._status_paths(helper._status_entries(target)) == set(
        manifest["allowlist"]
    )

    helper.execute(
        target,
        source_repo=source,
        manifest_path=manifest_path,
        verify_staged=True,
    )


def test_missing_target_is_rejected_before_any_write(pages_fixture):
    _, target, bundle, _, manifest = pages_fixture
    (target / "docs/source/insert.md").unlink()
    before = (target / "docs/source/page.md").read_text()
    with pytest.raises(helper.StageError, match="bounded-edit target is missing"):
        helper._plan_transformations(
            target, manifest, _sources(manifest, bundle), allow_section_drift=False
        )
    assert (target / "docs/source/page.md").read_text() == before


def test_missing_anchor_is_rejected(pages_fixture):
    _, target, bundle, _, manifest = pages_fixture
    _write(target / "docs/source/insert.md", "# Insert\n\nno anchor\n")
    with pytest.raises(helper.StageError, match="insertion anchor.*found 0"):
        helper._plan_transformations(
            target, manifest, _sources(manifest, bundle), allow_section_drift=False
        )


def test_duplicate_marker_is_rejected(pages_fixture):
    _, target, bundle, _, manifest = pages_fixture
    marker = f"<!-- BEGIN {helper.MARKER_PREFIX}: page -->"
    end = f"<!-- END {helper.MARKER_PREFIX}: page -->"
    _write(target / "docs/source/page.md", f"{marker}\n{end}\n{marker}\n{end}\n")
    with pytest.raises(helper.StageError, match="must be absent or appear once"):
        helper._plan_transformations(
            target, manifest, _sources(manifest, bundle), allow_section_drift=False
        )


def test_section_hash_drift_is_rejected(pages_fixture):
    _, target, bundle, _, manifest = pages_fixture
    _write(
        target / "docs/source/page.md",
        "# Page\n\n## Replace\ndrift\n\n## Next\nkeep\n",
    )
    with pytest.raises(helper.StageError, match="baseline section SHA256 drift"):
        helper._plan_transformations(
            target, manifest, _sources(manifest, bundle), allow_section_drift=False
        )


def test_planning_failure_leaves_all_targets_untouched(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    manifest = json.loads(manifest_path.read_text())
    manifest["targets"][1]["operations"][0]["anchor"] = "\n## Missing\n"
    manifest["targets"][1]["operations"][0]["anchor_sha256"] = _sha256(b"\n## Missing\n")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    before = (target / "docs/source/page.md").read_text()
    with pytest.raises(helper.StageError, match="insertion anchor.*found 0"):
        helper.execute(target, source_repo=source, manifest_path=manifest_path)
    assert (target / "docs/source/page.md").read_text() == before


def test_contradiction_search_rejects_stale_guidance(pages_fixture):
    _, target, bundle, _, manifest = pages_fixture
    plan = helper._plan_transformations(
        target, manifest, _sources(manifest, bundle), allow_section_drift=False
    )
    plan["docs/source/page.md"] += b"stale-guidance\n"
    with pytest.raises(helper.StageError, match="contradiction search failed"):
        helper._contradiction_report(target, plan, manifest)


def test_verify_staged_rejects_non_allowlisted_dirtiness(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    helper.execute(target, source_repo=source, manifest_path=manifest_path)
    _write(target / "docs/source/rogue.md", "not allowed\n")
    with pytest.raises(helper.StageError, match="exact allowlist"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            verify_staged=True,
        )


def test_verify_staged_rejects_allowlisted_surrounding_tampering(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    helper.execute(target, source_repo=source, manifest_path=manifest_path)
    page = target / "docs/source/page.md"
    _write(page, "tampered surroundings\n" + page.read_text())

    with pytest.raises(helper.StageError, match="differs from frozen HEAD plan"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            verify_staged=True,
        )


def test_verify_staged_rejects_index_tampering(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    helper.execute(target, source_repo=source, manifest_path=manifest_path)
    _run(target, "git", "add", "docs/source/page.md")

    with pytest.raises(helper.StageError, match="target index differs from HEAD"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            verify_staged=True,
        )


def test_verify_staged_rejects_intent_to_add(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    helper.execute(target, source_repo=source, manifest_path=manifest_path)
    _run(target, "git", "add", "-N", "docs/source/add.md")

    with pytest.raises(helper.StageError, match="intent-to-add"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            verify_staged=True,
        )


@pytest.mark.parametrize("flag", ["--assume-unchanged", "--skip-worktree"])
def test_strict_stage_rejects_hidden_index_flags(pages_fixture, flag):
    source, target, _, manifest_path, _ = pages_fixture
    unrelated = target / "docs/source/unrelated.md"
    _run(target, "git", "update-index", flag, unrelated.relative_to(target).as_posix())
    _write(unrelated, "hidden mutation\n")

    with pytest.raises(helper.StageError, match="noncanonical flag tag"):
        helper.execute(target, source_repo=source, manifest_path=manifest_path)


def test_verify_staged_rejects_allowlisted_symlink(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    helper.execute(target, source_repo=source, manifest_path=manifest_path)
    page = target / "docs/source/page.md"
    external = target.parent / "matching-page.md"
    external.write_bytes(page.read_bytes())
    page.unlink()
    page.symlink_to(external)

    with pytest.raises(helper.StageError, match="must be a regular file"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            verify_staged=True,
        )


def test_verify_staged_rejects_allowlisted_mode_tampering(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    helper.execute(target, source_repo=source, manifest_path=manifest_path)
    page = target / "docs/source/page.md"
    page.chmod(0o755)

    with pytest.raises(helper.StageError, match="mode drift"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            verify_staged=True,
        )


def test_reviewed_drift_packet_is_external_and_write_free(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    _write(
        target / "docs/source/page.md",
        "# Page\n\n## Replace\ndrift\n\n## Next\nkeep\n",
    )
    before = helper._target_fingerprint(target)
    packet = target.parent / "packet"

    helper.execute(
        target,
        source_repo=source,
        manifest_path=manifest_path,
        reviewed_drift=packet,
    )

    assert helper._target_fingerprint(target) == before
    assert (packet / "base/docs/source/page.md").is_file()
    assert (packet / "current/docs/source/page.md").is_file()
    assert (packet / "proposed/docs/source/page.md").is_file()
    assert (packet / "report.json").is_file()


def test_reviewed_drift_packet_rejects_other_registered_worktree(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    packet = source / "packet"

    with pytest.raises(helper.StageError, match="outside every registered worktree"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            reviewed_drift=packet,
        )
    assert not packet.exists()


def test_reviewed_drift_detects_ignored_target_write(pages_fixture, monkeypatch):
    source, target, _, manifest_path, _ = pages_fixture
    ignored = target / "docs/source/ignored.rcp07"
    _mutate_once_during_packet(monkeypatch, ignored, "ignored mutation\n")

    with pytest.raises(helper.StageError, match="modified the target worktree"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            reviewed_drift=target.parent / "packet",
        )
    assert _run(target, "git", "check-ignore", ignored.relative_to(target).as_posix())


def test_reviewed_drift_detects_dirty_non_allowlisted_mutation(
    pages_fixture, monkeypatch
):
    source, target, _, manifest_path, _ = pages_fixture
    unrelated = target / "docs/source/unrelated.md"
    _write(unrelated, "dirty before packet\n")
    _mutate_once_during_packet(monkeypatch, unrelated, "dirty during packet\n")

    with pytest.raises(helper.StageError, match="modified the target worktree"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            reviewed_drift=target.parent / "packet",
        )


def test_reviewed_drift_detects_target_root_metadata_write(pages_fixture, monkeypatch):
    source, target, _, manifest_path, _ = pages_fixture
    original = helper._write_packet_file
    original_mode = stat.S_IMODE(target.stat().st_mode)
    mutated = False

    def write_packet_file(*args, **kwargs):
        nonlocal mutated
        if not mutated:
            mutated = True
            target.chmod(original_mode ^ stat.S_IXOTH)
        return original(*args, **kwargs)

    monkeypatch.setattr(helper, "_write_packet_file", write_packet_file)
    with pytest.raises(helper.StageError, match="modified the target worktree"):
        helper.execute(
            target,
            source_repo=source,
            manifest_path=manifest_path,
            reviewed_drift=target.parent / "packet",
        )


def test_symbolic_branch_target_is_rejected(pages_fixture):
    source, target, _, manifest_path, _ = pages_fixture
    _run(target, "git", "switch", "-c", "pages-review")
    with pytest.raises(helper.StageError, match="target must be detached"):
        helper.execute(target, source_repo=source, manifest_path=manifest_path)


def test_public_manifest_matches_frozen_rcp07_inventory():
    manifest = helper._load_manifest(helper.MANIFEST_PATH)
    helper._validate_manifest(manifest, helper.BUNDLE_DIR)
    assert manifest["expected_baseline"] == "4833aa9341e19861297e330ff02aabfd8001935c"
    assert set(manifest["allowlist"]) == {
        "docs/source/configuration.md",
        "docs/source/examples/index.md",
        "docs/source/examples/io_outputs_and_sharding.md",
        "docs/source/modules/index.md",
        "docs/source/modules/outputs.md",
        "docs/source/running.md",
        "docs/source/tools/visualization.md",
        "docs/source/reference/input_parameters.md",
        "docs/source/reference/file_reference.md",
    }
    assert len(manifest["protected_blobs"]) == 4
    helper._load_sources(manifest, helper.BUNDLE_DIR)


def test_public_modules_fragments_preserve_table_rows_and_refresh_pdf_entry():
    manifest = helper._load_manifest(helper.MANIFEST_PATH)
    sources = helper._load_sources(manifest, helper.BUNDLE_DIR)
    item = next(
        item for item in manifest["targets"]
        if item["path"] == "docs/source/modules/index.md"
    )
    baseline = _run(
        REPO_ROOT, "git", "show", "origin/gh-pages:docs/source/modules/index.md"
    ) + "\n"
    result = helper._transform_target(
        item, baseline.encode(), sources, allow_section_drift=False
    ).decode()
    begin = f"<!-- BEGIN {helper.MARKER_PREFIX}: modules-index-output-count -->"
    marker = f"<!-- END {helper.MARKER_PREFIX}: modules-index-output-count -->"
    assert result.index(begin) < result.index("## Support Systems")
    assert result.index("| **Problem Generators** |") < result.index(marker)
    assert marker + "\n\n## Development Records" in result
    assert "[outputs.md](outputs.md)" in result
    assert "[boundaries.md](boundaries.md)" in result

    item = next(
        item for item in manifest["targets"]
        if item["path"] == "docs/source/modules/outputs.md"
    )
    baseline = _run(
        REPO_ROOT, "git", "show", "origin/gh-pages:docs/source/modules/outputs.md"
    ) + "\n"
    result = helper._transform_target(
        item, baseline.encode(), sources, allow_section_drift=False
    ).decode()
    begin = (
        f"<!-- BEGIN {helper.MARKER_PREFIX}: "
        "modules-outputs-implementation-entries -->"
    )
    marker = (
        f"<!-- END {helper.MARKER_PREFIX}: "
        "modules-outputs-implementation-entries -->"
    )
    assert result.index(begin) < result.index("## Implementation Entry Points")
    assert result.index("| `src/outputs/pdf.cpp` |") < result.index(marker)
    assert marker + "\n\n## See Also" in result
    assert "One- through four-dimensional histogram implementation" in result
    assert "One-/two-dimensional histogram implementation" not in result


def test_public_visualization_fragment_uses_header_path_and_scopes_header_api():
    fragment = (
        helper.BUNDLE_DIR
        / "fragments/docs/source/tools/visualization.io_outputs.md"
    ).read_text()
    assert (
        'pdf_header_path = "run/pdf_rho/node_00000000/simulation.header.pdf"'
        in fragment
    )
    assert "read_pdf.read_pdf_header(pdf_header_path, limits=limits)" in fragment
    assert "one header's intrinsic\nmetadata and declared shard identity" in fragment
    assert "Full readers discover siblings" in fragment


def test_public_file_reference_scopes_sphslice_header_api():
    fragment = (
        helper.BUNDLE_DIR
        / "fragments/docs/source/reference/file_reference.io_outputs.md"
    ).read_text()
    assert "selected shard's declared identity" in fragment
    assert "full\nreader discovers siblings" in fragment
