#!/usr/bin/env python3
"""Static contract tests for the bounded Q-032 reduced static-neutral carrier."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
PARTICLES = REPO_ROOT / "src/particles/particles.cpp"
SOURCE = REPO_ROOT / "src/pgen/tests/q032_reduced_static_neutral_local.cpp"
DECK = REPO_ROOT / "inputs/tests/pic_q032_reduced_static_neutral_runtime_local.athinput"
SCRIPT = (
    REPO_ROOT / "tst/scripts/particles/pic_q032_reduced_static_neutral_runtime_local.py"
)
SIDECAR = (
    REPO_ROOT / "tst/publication/readiness/"
    "q032_reduced_static_neutral_runtime_local_successor_2026-06-01.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_deck(path: Path) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            current = blocks.setdefault(line[1:-1], {})
            continue
        if current is not None and "=" in line:
            key, value = line.split("=", 1)
            current[key.strip()] = value.strip()
    return blocks


class Q032ReducedStaticNeutralRuntimeLocalTests(unittest.TestCase):
    def test_registration_is_additive_and_does_not_publish_blocked_open_name(self) -> None:
        pgen = (REPO_ROOT / "src/pgen/pgen.cpp").read_text(encoding="utf-8")
        self.assertEqual(pgen.count('compare("q032_reduced_static_neutral_local")'), 2)
        self.assertEqual(pgen.count("Q032ReducedStaticNeutralLocal(pin, false);"), 1)
        self.assertEqual(pgen.count("Q032ReducedStaticNeutralLocal(pin, true);"), 1)
        self.assertNotIn('compare("q032_plotnikov_damped_crsi_open")', pgen)
        self.assertIn(
            "void Q032ReducedStaticNeutralLocal(ParameterInput *pin, const bool restart);",
            (REPO_ROOT / "src/pgen/pgen.hpp").read_text(encoding="utf-8"),
        )
        self.assertIn(
            "pgen/tests/q032_reduced_static_neutral_local.cpp",
            (REPO_ROOT / "src/CMakeLists.txt").read_text(encoding="utf-8"),
        )
        self.assertIn(
            "`q032_reduced_static_neutral_local`",
            (REPO_ROOT / "src/pgen/AGENTS.md").read_text(encoding="utf-8"),
        )

    def test_deck_is_exact_thin_isolated_mechanics_carrier(self) -> None:
        blocks = _parse_deck(DECK)
        self.assertEqual(
            {key: blocks["mesh"][key] for key in ("nx1", "nx2", "nx3")},
            {"nx1": "32", "nx2": "4", "nx3": "1"},
        )
        self.assertEqual(
            {key: blocks["meshblock"][key] for key in ("nx1", "nx2", "nx3")},
            {"nx1": "32", "nx2": "4", "nx3": "1"},
        )
        for axis in ("ix1_bc", "ox1_bc", "ix2_bc", "ox2_bc"):
            self.assertEqual(blocks["mesh"][axis], "periodic")
        self.assertEqual(blocks["mesh_refinement"]["refinement"], "none")
        particles = blocks["particles"]
        self.assertEqual(particles["pic_physical_mode"], "extended_mhd_pic")
        self.assertEqual(particles["pic_wave_damping_mode"], "ion_neutral_friction")
        self.assertEqual(particles["pic_ion_neutral_collision_rate"], "0.7")
        self.assertEqual(particles["pic_feedback_mode"], "test_particle")
        self.assertEqual(particles["pic_enable_2d3v"], "true")
        self.assertEqual(particles["ppc"], "0.0")
        self.assertEqual(particles["deposit_moments"], "false")
        self.assertEqual(particles["couple_moments_to_mhd"], "false")
        self.assertEqual(blocks["problem"]["pgen_name"],
                         "q032_reduced_static_neutral_local")

    def test_generator_rejects_broadened_carriers_and_delegates_mhd_initialization(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        self.assertIn("q032_reduced_static_neutral_local rejects AMR/SMR", source)
        self.assertIn("requires periodic active ", source)
        self.assertIn("32x4x1 2D3V x1-parallel carrier with 32x4x1 serial or ", source)
        self.assertIn("16x4x1 two-way x1 meshblocks", source)
        self.assertIn('Q032RequireString(pin, "particles", "pic_feedback_mode", '
                      '"test_particle");', source)
        self.assertIn("requires the exact Newtonian ", source)
        self.assertIn("single-fluid MHD source task path", source)
        self.assertIn("pmbp->pcoord->is_special_relativistic", source)
        self.assertIn("pmbp->pcoord->is_general_relativistic", source)
        self.assertIn("pmbp->pcoord->is_dynamical_relativistic", source)
        self.assertIn('pin->DoesBlockExist("radiation")', source)
        self.assertIn('pin->DoesBlockExist("ion-neutral")', source)
        self.assertIn('Q032RequireString(pin, block, "plotnikov_qualification", '
                      '"not_claimed");', source)
        self.assertIn("LinearWave(pin, restart);", source)

    def test_shared_parser_rejects_non_newtonian_or_alternate_damping_paths(self) -> None:
        particles = PARTICLES.read_text(encoding="utf-8")
        self.assertIn("the Newtonian single-fluid MHD source task path", particles)
        for reason in (
            "<radiation> task paths",
            "<ion-neutral> alternate task lists",
            "<hydro> compositions",
            "dynamical-GR <adm>/<z4c> task paths",
            "<coord>/special_rel=true",
            "<coord>/general_rel=true",
            "dynamical-relativistic coordinates",
        ):
            self.assertIn(reason, particles)

    def test_sidecar_binds_only_owned_artifacts_and_retains_nonclaims(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-032")
        self.assertEqual(sidecar["qualification_effect"], "none")
        self.assertFalse(sidecar["claim_closure"])
        self.assertEqual(sidecar["frontier_authorization"], "not_bound")
        self.assertEqual(
            sidecar["generator"]["carrier"]["fluid_task_path"],
            "newtonian_single_fluid_mhd_source_task_only",
        )
        self.assertEqual(
            sidecar["generator"]["carrier"]["mesh"],
            {"nx1": 32, "nx2": 4, "nx3": 1},
        )
        self.assertEqual(
            sidecar["generator"]["carrier"]["meshblock"],
            {"nx1": "32 serial or 16 two-way x1", "nx2": 4, "nx3": 1},
        )
        expected_paths = {
            str(PARTICLES.relative_to(REPO_ROOT)),
            str(SOURCE.relative_to(REPO_ROOT)),
            str(DECK.relative_to(REPO_ROOT)),
            str(SCRIPT.relative_to(REPO_ROOT)),
            str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        }
        self.assertEqual(set(sidecar["owned_artifact_bindings"]), expected_paths)
        for relative, expected_sha256 in sidecar["owned_artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected_sha256)
        self.assertTrue(
            any("Plotnikov" in item for item in sidecar["explicitly_not_claimed"])
        )
        self.assertTrue(
            any("Frontier" in item for item in sidecar["explicitly_not_claimed"])
        )
        self.assertTrue(
            any(
                "alternate-task-path damping support" in item
                for item in sidecar["explicitly_not_claimed"]
            )
        )
        for rejection in (
            "special-relativistic MHD rejected",
            "general-relativistic MHD rejected",
            "radiation alternate task path rejected",
            "ion-neutral alternate task list rejected",
            "dynamical-GR alternate task path rejected",
        ):
            self.assertIn(
                rejection,
                sidecar["source_local_runtime_probe"]["negative_contract"],
            )
        runtime_result = sidecar["source_local_runtime_probe"]["result"]
        self.assertEqual(runtime_result["status"], "pass")
        self.assertEqual(runtime_result["qualification_effect"], "none")
        self.assertEqual(runtime_result["plotnikov_qualification"], "not_claimed")
        self.assertLessEqual(max(runtime_result["relative_errors"].values()), 1.0e-7)
        self.assertEqual(len(runtime_result["executable_sha256"]), 64)

    def test_runtime_probe_includes_mechanics_comparison_and_rejection_paths(self) -> None:
        script = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("'problem/pgen_name=linear_wave'", script)
        self.assertIn("'particles/pic_wave_damping_mode=off'", script)
        self.assertIn("'mesh_refinement/refinement=static'", script)
        self.assertIn("'mesh/ix1_bc=outflow'", script)
        self.assertIn("q032_reduced_static_neutral_local/deck_role=plotnikov_evidence",
                      script)
        for guard in (
            "rejects_wrapper_special_relativistic",
            "rejects_special_relativistic",
            "rejects_general_relativistic",
            "rejects_radiation_task_path",
            "rejects_ion_neutral_task_path",
            "rejects_dynamical_gr_task_path",
        ):
            self.assertIn(guard, script)
        self.assertIn("_require_equal_snapshots('initial control/damped'", script)
        self.assertIn("'final control/damped snapshot times differ'", script)
        self.assertIn("'qualification_effect': 'none'", script)
        self.assertIn("'plotnikov_qualification': 'not_claimed'", script)
        self.assertIn("ATHENA_Q032_NPROC", script)
        self.assertIn("ATHENA_Q032_LAUNCHER", script)
        self.assertIn("['meshblock/nx1=16']", script)
        self.assertIn("'configured_ranks': _rank_count()", script)


if __name__ == "__main__":
    unittest.main()
