#!/usr/bin/env python3
"""Regression tests for PIC engineering-artifact taxonomy and safety rails."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[2]
PUBLICATION_DIR = REPO_ROOT / "tst" / "publication"
sys.path.insert(0, str(PUBLICATION_DIR))

import pic_publication_manifest as manifest  # noqa: E402
import run_pic_publication_suite as suite  # noqa: E402
import run_pic_shock_scan as shock_scan  # noqa: E402
from artifact_lineage import write_companion_manifest  # noqa: E402


class PicArtifactTaxonomyTests(unittest.TestCase):
    def test_manifest_cases_are_explicit_engineering_proxies(self) -> None:
        legacy_groups = set(manifest.LEGACY_GROUP_ALIASES)
        for case in manifest.ARTIFACT_CASES:
            self.assertEqual(case.evidence_class, "engineering_proxy")
            self.assertNotEqual(case.physical_model, "exploratory_unqualified")
            self.assertNotIn(case.group, legacy_groups)

    def test_legacy_group_aliases_select_the_same_cases(self) -> None:
        for legacy, canonical in manifest.LEGACY_GROUP_ALIASES.items():
            old_ids = {
                case.case_id for case in manifest.select_cases([legacy], "all")
            }
            new_ids = {
                case.case_id for case in manifest.select_cases([canonical], "all")
            }
            self.assertEqual(old_ids, new_ids)

    def test_visible_crsi_proxy_labels_do_not_claim_delta_f(self) -> None:
        visible_sources = (
            PUBLICATION_DIR / "pic_publication_manifest.py",
            PUBLICATION_DIR / "plot_pic_publication_figures.py",
        )
        for path in visible_sources:
            text = path.read_text()
            self.assertNotIn("delta-f noise suppression", text)
            self.assertNotIn("CRSI delta-f polarization growth", text)

    def test_archive_outputs_are_checksummed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            case_dir = Path(tmp) / "case"
            raw = repo / "raw" / "value.bin"
            raw.parent.mkdir(parents=True)
            raw.write_bytes(b"artifact-lineage")
            records = suite._archive_outputs(repo, case_dir, ("raw/*.bin",))
            self.assertEqual(len(records), 1)
            self.assertEqual(
                records[0]["sha256"],
                hashlib.sha256(b"artifact-lineage").hexdigest(),
            )
            self.assertTrue((case_dir / records[0]["archive_path"]).is_file())

    def test_scan_defaults_and_hpc_safety(self) -> None:
        self.assertEqual(
            shock_scan._resolve_default_lists("local"),
            ("0.8,1.0,1.2", "2,4,8", "32,48,64"),
        )
        self.assertEqual(
            shock_scan._resolve_default_lists("hpc"),
            ("1.0,1.2,1.4", "8,16,24", "256,384,512"),
        )
        with self.assertRaises(ValueError):
            shock_scan._resolve_deck(REPO_ROOT, "local", "hpc")

        proc = subprocess.run(
            [
                sys.executable,
                str(PUBLICATION_DIR / "run_pic_shock_scan.py"),
                "--tier",
                "hpc",
                "--run",
            ],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("--tier hpc is emit-only", proc.stderr)

        with tempfile.TemporaryDirectory() as tmp:
            emit_proc = subprocess.run(
                [
                    sys.executable,
                    str(PUBLICATION_DIR / "run_pic_shock_scan.py"),
                    "--tier",
                    "hpc",
                    "--max-cases",
                    "1",
                    "--output-dir",
                    tmp,
                    "--run-id",
                    "emit_only",
                ],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(emit_proc.returncode, 0, emit_proc.stderr)
            script = (
                Path(tmp)
                / "emit_only"
                / "reviews"
                / "shock_step4"
                / "scan_commands.sh"
            ).read_text()
            self.assertIn("planning material only", script)
            self.assertIn("exit 2", script)
            self.assertIn("# mpiexec", script)
            lineage = json.loads(
                (
                    Path(tmp)
                    / "emit_only"
                    / "reviews"
                    / "shock_step4"
                    / "scan_artifact_manifest.json"
                ).read_text()
            )
            self.assertTrue(lineage["generator_source_sha256"])
            self.assertTrue(lineage["git_commit"])
            self.assertEqual(len(lineage["outputs"]), 3)

        with tempfile.TemporaryDirectory() as tmp:
            local_script = Path(tmp) / "scan_commands.sh"
            shock_scan._write_command_script(local_script, REPO_ROOT, [], "local")
            self.assertIn("planning material only", local_script.read_text())
            self.assertIn("exit 2", local_script.read_text())
            self.assertEqual(local_script.stat().st_mode & 0o111, 0)

    def test_shock_scaffold_direct_run_requires_local_acknowledgement(self) -> None:
        proc = subprocess.run(
            [
                sys.executable,
                str(PUBLICATION_DIR / "run_pic_parallel_shock_benchmark.py"),
                "--run",
            ],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("--run requires --local-exploratory-run", proc.stderr)

        for variable in ("SLURM_JOB_ID", "SLURM_JOBID"):
            env = dict(os.environ)
            env[variable] = "taxonomy-test"
            guarded = subprocess.run(
                [
                    sys.executable,
                    str(PUBLICATION_DIR / "run_pic_parallel_shock_benchmark.py"),
                    "--run",
                    "--local-exploratory-run",
                ],
                cwd=str(REPO_ROOT),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(guarded.returncode, 0)
            self.assertIn("direct execution is local-only", guarded.stderr)

    def test_scan_direct_run_requires_local_acknowledgement(self) -> None:
        proc = subprocess.run(
            [
                sys.executable,
                str(PUBLICATION_DIR / "run_pic_shock_scan.py"),
                "--tier",
                "local",
                "--run",
            ],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("local execution requires --local-exploratory-run", proc.stderr)

    def test_suite_direct_run_requires_local_acknowledgement(self) -> None:
        proc = subprocess.run(
            [
                sys.executable,
                str(PUBLICATION_DIR / "run_pic_publication_suite.py"),
                "--skip-build",
                "--tier",
                "local",
                "--cases",
                "bell_growth_publication",
            ],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("local execution requires --local-exploratory-run", proc.stderr)

    def test_suite_rejects_hpc_and_slurm_direct_execution(self) -> None:
        base = [
            sys.executable,
            str(PUBLICATION_DIR / "run_pic_publication_suite.py"),
            "--skip-build",
            "--cases",
            "bell_growth_publication",
            "--local-exploratory-run",
        ]
        for tier in ("hpc", "all"):
            proc = subprocess.run(
                base + ["--tier", tier],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("direct HPC/all-tier execution is prohibited", proc.stderr)

        for variable in ("SLURM_JOB_ID", "SLURM_JOBID"):
            env = dict(os.environ)
            env[variable] = "taxonomy-test"
            proc = subprocess.run(
                base + ["--tier", "local"],
                cwd=str(REPO_ROOT),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("direct execution is local-only", proc.stderr)

    def test_scan_rejects_slurm_direct_execution(self) -> None:
        for variable in ("SLURM_JOB_ID", "SLURM_JOBID"):
            env = dict(os.environ)
            env[variable] = "taxonomy-test"
            proc = subprocess.run(
                [
                    sys.executable,
                    str(PUBLICATION_DIR / "run_pic_shock_scan.py"),
                    "--tier",
                    "local",
                    "--run",
                    "--local-exploratory-run",
                    "--max-cases",
                    "1",
                ],
                cwd=str(REPO_ROOT),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("direct execution is local-only", proc.stderr)

    def test_reproduction_plot_bundle_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(
                [
                    sys.executable,
                    str(PUBLICATION_DIR / "plot_pic_publication_figures.py"),
                    "--run-dir",
                    tmp,
                    "--bundle",
                    "sun_bai_2023_reproduction",
                ],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("No checked-in Sun and Bai", proc.stderr)

    def test_dry_run_summary_and_manifest_record_taxonomy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(
                [
                    sys.executable,
                    str(PUBLICATION_DIR / "run_pic_publication_suite.py"),
                    "--dry-run",
                    "--skip-build",
                    "--tier",
                    "local",
                    "--cases",
                    "bell_growth_publication",
                    "--run-root",
                    tmp,
                    "--run-id",
                    "taxonomy",
                ],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            run_dir = Path(tmp) / "taxonomy"
            summary = json.loads((run_dir / "summary.json").read_text())
            artifact_manifest = json.loads(
                (run_dir / "artifact_manifest.json").read_text()
            )
            lineage_manifest = json.loads((run_dir / "lineage_manifest.json").read_text())
            self.assertEqual(summary["evidence_class"], "engineering_proxy")
            self.assertTrue(summary["not_sun_bai_reproduction"])
            self.assertTrue(summary["cases"][0]["module_source_sha256"])
            self.assertEqual(
                artifact_manifest["evidence_class"], "engineering_proxy"
            )
            self.assertTrue(artifact_manifest["files"])
            self.assertEqual(lineage_manifest["schema_version"], 1)
            self.assertTrue(lineage_manifest["generator_source_sha256"])
            self.assertTrue(lineage_manifest["lineage_helper_source_sha256"])

    def test_dry_run_does_not_delete_matching_existing_artifact(self) -> None:
        case = next(
            item
            for item in manifest.ARTIFACT_CASES
            if item.case_id == "bell_growth_publication"
        )
        seed = REPO_ROOT / case.output_globs[0].replace("*", "taxonomy")
        seed.parent.mkdir(parents=True, exist_ok=True)
        seed.write_bytes(b"dry-run-must-not-delete")
        try:
            with tempfile.TemporaryDirectory() as tmp:
                proc = subprocess.run(
                    [
                        sys.executable,
                        str(PUBLICATION_DIR / "run_pic_publication_suite.py"),
                        "--dry-run",
                        "--skip-build",
                        "--tier",
                        "local",
                        "--cases",
                        case.case_id,
                        "--run-root",
                        tmp,
                        "--run-id",
                        "dry_run_preserves_artifacts",
                    ],
                    cwd=str(REPO_ROOT),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(proc.returncode, 0, proc.stderr)
                self.assertTrue(seed.is_file())
        finally:
            seed.unlink(missing_ok=True)

    def test_retired_runbook_is_non_operational(self) -> None:
        text = (PUBLICATION_DIR / "PIC_LARGE_MACHINE_VALIDATION.md").read_text()
        self.assertIn("Historical evidence only", text)
        self.assertNotIn("```", text)
        self.assertNotIn("Treat the branch as large-machine validated", text)

    def test_root_legacy_pic_documents_are_historical_only(self) -> None:
        for path in sorted(REPO_ROOT.glob("AGENT_PIC*.md")):
            first_lines = "\n".join(path.read_text().splitlines()[:4])
            self.assertIn("HISTORICAL ONLY", first_lines, str(path))
            self.assertIn("PIC_PRODUCTION_READINESS_PLAN.md", first_lines, str(path))

    def test_nested_orientation_docs_use_scaffold_taxonomy(self) -> None:
        inputs_text = (REPO_ROOT / "inputs" / "AGENTS.md").read_text()
        pgen_text = (REPO_ROOT / "src" / "pgen" / "AGENTS.md").read_text()
        self.assertNotIn("Heavy production shock deck", inputs_text)
        self.assertNotIn("Section 5.4-aligned parallel-shock benchmark", inputs_text)
        self.assertNotIn("physics-benchmark pgen", pgen_text)
        self.assertIn("unqualified engineering scaffolds", pgen_text)

    def test_standalone_helpers_write_checksummed_lineage_companions(self) -> None:
        retained_helpers = (
            "analyze_pic_shock_storyline.py",
            "plot_pic_shock_storyline.py",
            "plot_pic_parallel_shock_section54.py",
            "run_pic_parallel_shock_benchmark.py",
            "run_pic_shock_scan.py",
            "plot_pic_shock_paper_style.py",
            "plot_pic_publication_figures.py",
            "compute_pic_publication_metrics.py",
            "check_pic_parallel_shock_feedback_sensitivity.py",
        )
        for helper in retained_helpers:
            text = (PUBLICATION_DIR / helper).read_text()
            self.assertIn("write_companion_manifest(", text, helper)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.bin"
            output = root / "figure.png"
            source.write_bytes(b"source")
            output.write_bytes(b"output")
            executable = root / "athena"
            executable.write_bytes(b"executable")
            companion = root / "artifact_manifest.json"
            write_companion_manifest(
                companion,
                generator="test_pic_artifact_taxonomy.py",
                physical_model="test-fixture",
                inputs=(source,),
                outputs=(output,),
                executables=(executable,),
            )
            data = json.loads(companion.read_text())
            self.assertEqual(data["schema_version"], 1)
            self.assertTrue(data["generator_source_sha256"])
            self.assertTrue(data["lineage_helper_source_sha256"])
            self.assertTrue(data["git_commit"])
            self.assertEqual(
                data["inputs"][0]["sha256"], hashlib.sha256(b"source").hexdigest()
            )
            self.assertEqual(
                data["outputs"][0]["sha256"], hashlib.sha256(b"output").hexdigest()
            )
            self.assertEqual(
                data["executables"][0]["sha256"],
                hashlib.sha256(b"executable").hexdigest(),
            )
            with self.assertRaises(FileNotFoundError):
                write_companion_manifest(
                    root / "missing_input_manifest.json",
                    generator="test_pic_artifact_taxonomy.py",
                    physical_model="test-fixture",
                    inputs=(root / "missing.bin",),
                    outputs=(output,),
                )
            with self.assertRaises(FileNotFoundError):
                write_companion_manifest(
                    root / "missing_generator_manifest.json",
                    generator="missing_generator.py",
                    physical_model="test-fixture",
                    inputs=(source,),
                    outputs=(output,),
                )
            with mock.patch("artifact_lineage.subprocess.run") as git_run:
                git_run.return_value.returncode = 1
                git_run.return_value.stdout = ""
                with self.assertRaises(RuntimeError):
                    write_companion_manifest(
                        root / "missing_git_manifest.json",
                        generator="test_pic_artifact_taxonomy.py",
                        physical_model="test-fixture",
                        inputs=(source,),
                        outputs=(output,),
                    )

    def test_lineage_clean_worktree_marker_is_reachable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            repo = root / "repo"
            repo.mkdir()
            tracked = repo / "tracked.txt"
            tracked.write_text("tracked\n")
            subprocess.run(
                ["git", "init", "-b", "main"],
                cwd=str(repo),
                capture_output=True,
                text=True,
                check=True,
            )
            subprocess.run(
                ["git", "add", "tracked.txt"],
                cwd=str(repo),
                capture_output=True,
                text=True,
                check=True,
            )
            subprocess.run(
                [
                    "git",
                    "-c",
                    "user.name=PIC Taxonomy Test",
                    "-c",
                    "user.email=pic-taxonomy@example.invalid",
                    "commit",
                    "-m",
                    "Initialize clean fixture",
                ],
                cwd=str(repo),
                capture_output=True,
                text=True,
                check=True,
            )
            source = root / "source.bin"
            output = root / "output.bin"
            source.write_bytes(b"source")
            output.write_bytes(b"output")
            companion = root / "clean_manifest.json"
            with mock.patch("artifact_lineage.REPO_ROOT", repo):
                write_companion_manifest(
                    companion,
                    generator="test_pic_artifact_taxonomy.py",
                    physical_model="test-fixture",
                    inputs=(source,),
                    outputs=(output,),
                )
            data = json.loads(companion.read_text())
            self.assertEqual(data["git_branch"], "main")
            self.assertEqual(data["git_status"], ["<clean>"])

    def test_lineage_callers_do_not_prefilter_required_files(self) -> None:
        storyline = (PUBLICATION_DIR / "plot_pic_shock_storyline.py").read_text()
        figures = (PUBLICATION_DIR / "plot_pic_publication_figures.py").read_text()
        scan = (PUBLICATION_DIR / "run_pic_shock_scan.py").read_text()
        self.assertIn("inputs.append(Path(args.scan_summary).resolve())", storyline)
        self.assertNotIn(
            "[path for path in generated_files if path.is_file()]", figures
        )
        self.assertIn("executables=(athena,) if args.run else ()", scan)

    def test_legacy_decks_are_visibly_labeled(self) -> None:
        decks = sorted((REPO_ROOT / "inputs" / "tests").glob("*publication*.athinput"))
        decks.append(
            REPO_ROOT / "inputs" / "tests" / "pic_parallel_shock_section54_hpc.athinput"
        )
        decks.extend(
            [
                REPO_ROOT
                / "inputs"
                / "tests"
                / "pic_parallel_shock_coarse_uniform.athinput",
                REPO_ROOT
                / "inputs"
                / "tests"
                / "pic_parallel_shock_fine_uniform.athinput",
                REPO_ROOT
                / "inputs"
                / "tests"
                / "pic_parallel_shock_amr_fiducial.athinput",
            ]
        )
        for deck in decks:
            first_line = deck.read_text().splitlines()[0]
            self.assertTrue(
                first_line.startswith(("# LEGACY NAME:", "# UNQUALIFIED ENGINEERING")),
                str(deck),
            )


if __name__ == "__main__":
    unittest.main()
