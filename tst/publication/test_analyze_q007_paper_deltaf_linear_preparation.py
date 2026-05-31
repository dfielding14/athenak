#!/usr/bin/env python3
"""Focused tests for bounded Q-007 true-delta-f source-local preparation."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import stat
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from tst.publication import analyze_q007_paper_deltaf_linear_preparation as q007
from tst.publication.pvtk_particles import ParticleVTKData


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q007_paper_deltaf_linear_source_local_preparation_successor_v2_2026-05-31.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Sparse values independently pinned from the source-local loading equations.
# Keep these literal: deriving them through q007 would mirror analyzer mistakes.
_PINNED_STARTUP_LOADING_SAMPLES = (
    (
        "crsi", 0, 0, (5.0, 0.5, 0.0),
        (0.43816292801599277, 0.10437500304887135, -1.2245219489485344),
        (0.4381587841856521, 0.10437401594504127, -1.224510368299236),
        7.51801921412431e-10, 0.9999659541392419,
    ),
    (
        "crsi", 7, 7, (5.0, 0.5, 0.0),
        (23165.044439454363, 5518.156442269891, -64738.71601352069),
        (100.74657173369104, 23.998889589955287, -281.55368809457985),
        4.3305272340061035e-06, 3.897894082769363e-11,
    ),
    (
        "crsi", 8, 0, (15.0, 0.5, 0.0),
        (0.43816292801599277, 0.10437500304887135, -1.2245219489485344),
        (0.4381587841856521, 0.10437401594504127, -1.224510368299236),
        7.51801921412431e-10, 0.9999659541392419,
    ),
    (
        "crsi", 256, 0, (5.0, 1.5, 0.0),
        (0.43816292801599277, 0.10437500304887135, -1.2245219489485344),
        (0.4381587841856521, 0.10437401594504127, -1.224510368299236),
        7.51801921412431e-10, 0.9999659541392419,
    ),
    (
        "crsi", 1024, 0, (5.0, 0.5, 0.0),
        (-0.43816292801599277, -0.10437500304887135, 1.2245219489485344),
        (-0.4381587841856521, -0.10437401594504127, 1.224510368299236),
        7.51801921412431e-10, 0.9999659541392419,
    ),
    (
        "crsi", 2048, 0, (5.0, 0.5, 0.0),
        (-1.2201851496905332, 0.30091217853927954, -0.35062081724728283),
        (-1.220173610055568, 0.3009093327279715, -0.35061750132736424),
        7.51801921412431e-10, 0.9999659541392419,
    ),
    (
        "crpai_prolate", 0, 0, (10.0, 0.5, 0.0),
        (0.4243459833398164, 0.5201804359676256, 1.1325130683713631),
        (0.4243459829312082, 0.5201804354667372, 1.1325130672808517),
        8.767272487179971e-10, 0.9999702772946475,
    ),
    (
        "crpai_oblate", 0, 0, (10.0, 0.5, 0.0),
        (0.4243459833398164, 0.5098798332751975, 1.1100870670174747),
        (0.4243459829455659, 0.5098798328014794, 1.1100870659861173),
        8.767272487179971e-10, 0.9999702772946475,
    ),
)


def _startup_particle_payload(case: str) -> ParticleVTKData:
    expected = q007._expected_startup_loading(case)
    return ParticleVTKData(
        points=expected["points"].copy(),
        scalars={
            "ptag": expected["tags"].copy(),
            "species": expected["species"].copy(),
            "macro_weight": expected["weights"].copy(),
            "deltaf_f0": expected["f0"].copy(),
            "deltaf_weight": np.zeros(q007.SOURCE_LOCAL_PARTICLE_TOTAL),
        },
        vectors={"vel": expected["velocities"].astype(np.float32).astype(np.float64)},
    )


def _startup_mhd_payload(case: str) -> dict[str, object]:
    expected = q007._CASES[case]
    length = float(expected["x1max"])
    x1v = (np.arange(32) + 0.5) * (length / 32.0)
    state = np.array([
        q007.four_branch_alfven_state(
            x1, length, float(expected["paper_seed_amplitude"]), int(expected["wave_seed"])
        )
        for x1 in x1v
    ])
    shape = (1, 4, 32)
    return {
        "MaxLevel": 0,
        "Time": 0.0,
        "NumCycles": 0,
        "x1f": np.linspace(0.0, length, 33),
        "x1v": x1v,
        "x2f": np.arange(5, dtype=np.float64),
        "x2v": np.arange(4, dtype=np.float64) + 0.5,
        "x3f": np.arange(2, dtype=np.float64),
        "x3v": np.array([0.5]),
        "dens": np.broadcast_to(np.ones(32), shape).copy(),
        "velx": np.broadcast_to(
            np.full(32, float(expected["source_local_gas_vx"])), shape
        ).copy(),
        "vely": np.broadcast_to(state[:, 2], shape).copy(),
        "velz": np.broadcast_to(state[:, 3], shape).copy(),
        "bcc1": np.broadcast_to(np.ones(32), shape).copy(),
        "bcc2": np.broadcast_to(state[:, 0], shape).copy(),
        "bcc3": np.broadcast_to(state[:, 1], shape).copy(),
    }


def _write_runtime_invocation_tree(root: Path) -> dict[str, object]:
    executable = str(root.parent / "pin/bin/athena")
    pinned = {
        "pinned_executable_path": executable,
        "pinned_executable_root": str(root.parent / "pin"),
        "pinned_executable_root_inventory_sha256": "1" * 64,
        "pinned_executable_sha256": "2" * 64,
    }
    cases = {
        "crsi": ("crsi", "pic_q007_paper_crsi_linear_preparation.athinput", [], 0, range(4)),
        "crpai-prolate": (
            "crpai_prolate",
            "pic_q007_paper_crpai_linear_prolate_preparation.athinput",
            [],
            0,
            range(2),
        ),
        "crpai-oblate": (
            "crpai_oblate",
            "pic_q007_paper_crpai_linear_oblate_preparation.athinput",
            [],
            0,
            range(2),
        ),
        "negative-crpai-nlim1": (
            "crpai_prolate",
            "pic_q007_paper_crpai_linear_prolate_preparation.athinput",
            ["time/nlim=1"],
            1,
            range(0),
        ),
    }
    deck_root = root / "decks"
    deck_root.mkdir(parents=True)
    for case, name, _, _, _ in cases.values():
        (deck_root / name).write_bytes(q007.DECKS[case].read_bytes())
    for label, (case, deck_name, overrides, returncode, cycles) in cases.items():
        run = root / label
        run.mkdir()
        stdout = (
            "<time>/nlim does not match the Q-007 preparation contract\n"
            if label == "negative-crpai-nlim1"
            else "bounded\n"
        )
        (run / "athena.stdout.txt").write_text(stdout, encoding="utf-8")
        (run / "athena.stderr.txt").write_text("", encoding="utf-8")
        (run / "returncode.txt").write_text(f"{returncode}\n", encoding="utf-8")
        basename = q007._CASES[case]["basename"]
        for cycle in cycles:
            (run / "bin").mkdir(exist_ok=True)
            (run / "pvtk").mkdir(exist_ok=True)
            (run / "bin" / f"{basename}.mhd_w_bcc.{cycle:05d}.bin").write_text(
                "bounded\n", encoding="utf-8"
            )
            (run / "pvtk" / f"{basename}.prtcl_all.{cycle:05d}.part.vtk").write_text(
                "bounded\n", encoding="utf-8"
            )
        deck = deck_root / deck_name
        invocation = {
            "schema_version": 1,
            "argv": [executable, "-i", str(deck), *overrides],
            "cwd": str(run),
            "deck": {"path": str(deck), "sha256": _sha256(deck)},
            "executable_realpath": executable,
            "executable_sha256": pinned["pinned_executable_sha256"],
            "frontier_used": False,
            "kronos_used": False,
            "launch_style": "direct_serial_host_execution",
            "mpi_used": False,
            "overrides": overrides,
            "pinned_executable_root": pinned["pinned_executable_root"],
            "pinned_executable_root_inventory_sha256":
                pinned["pinned_executable_root_inventory_sha256"],
            "produced_file_sha256": {
                path.relative_to(run).as_posix(): _sha256(path)
                for path in sorted(run.rglob("*")) if path.is_file()
            },
            "returncode": returncode,
            "selected_environment": {"PYTHONDONTWRITEBYTECODE": "1"},
            "slurm_used": False,
        }
        (run / "invocation.json").write_text(
            json.dumps(invocation, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return pinned


class Q007PaperDeltaFLinearPreparationTests(unittest.TestCase):
    def test_crsi_static_paper_mapping(self) -> None:
        contract = q007.analytical_contract()
        crsi = contract["crsi"]
        self.assertAlmostEqual(crsi["k0"], 1.0 / 300.0)
        self.assertAlmostEqual(crsi["lambda0"], 600.0 * math.pi)
        self.assertGreater(crsi["q1_at_k0"], 0.0)
        self.assertGreater(crsi["q2_at_k0"], 0.0)
        self.assertGreater(crsi["forward_growth_at_k0"], 0.0)
        self.assertLess(crsi["backward_growth_at_minus_k0"], 0.0)
        self.assertTrue(contract["full_q1_q2_dispersion_oracle_prepared"])
        self.assertFalse(contract["runtime_growth_fit_claimed"])
        self.assertFalse(contract["qualifying_evidence"])

    def test_crpai_distribution_and_signed_branch_mapping(self) -> None:
        isotropic = q007.isotropic_kappa_distribution(1.0, 0.0, 300.0, 1.75)
        for role, xi, branch in (("prolate", 0.99, "-1"), ("oblate", 1.01, "1")):
            with self.subTest(role=role):
                anisotropic = q007.anisotropic_kappa_distribution(
                    1.0, 0.0, 0.0, 0.0, 300.0, 1.75, xi
                )
                self.assertAlmostEqual(anisotropic / isotropic, xi * xi)
                case = q007.analytical_contract()["crpai"][role]
                self.assertAlmostEqual(case["athenak_transverse_anisotropy_scale"],
                                       1.0 / xi)
                self.assertEqual(case["unstable_signed_branch"], branch)
                self.assertGreater(case["signed_branch_growth"][branch], 0.0)
                self.assertEqual(
                    case["handedness_mapping"],
                    "blocked_pending_manuscript_text_caption_review",
                )
        self.assertFalse(q007.analytical_contract()["crpai_handedness_claimed"])

    def test_three_decks_freeze_weighted_loading_and_bounded_replay(self) -> None:
        decks = q007.validate_decks()
        self.assertEqual(
            [item["case"] for item in decks],
            ["crsi", "crpai_prolate", "crpai_oblate"],
        )
        self.assertEqual(
            [item["cycle_zero_only"] for item in decks], [False, True, True]
        )
        self.assertTrue(all(item["true_deltaf_parser_path"] for item in decks))
        self.assertTrue(all(item["exact_isothermal_mhd"] for item in decks))
        self.assertTrue(
            all(item["source_local_particles_total"] == 262144 for item in decks)
        )
        self.assertTrue(
            all(item["paper_particles_per_cell_total"] == 2048 for item in decks)
        )
        self.assertTrue(all(item["physical_loading_implemented"] for item in decks))
        self.assertTrue(
            all(item["random_phase_wave_spectrum_implemented"] for item in decks)
        )
        self.assertEqual(
            [item["runtime_evolution_admitted"] for item in decks], [True, False, False]
        )
        self.assertTrue(all(not item["qualifying_evidence"] for item in decks))

    def test_deck_drift_fails_closed(self) -> None:
        mutations = (
            ("crsi", "nlim       = 2", "nlim       = 3", "time/nlim"),
            (
                "crsi",
                "couple_moments_energy_to_mhd      = false",
                "couple_moments_energy_to_mhd      = true",
                "couple_moments_energy_to_mhd",
            ),
            (
                "crsi",
                "couple_moments_momentum_coeff     = 1.0",
                "couple_moments_momentum_coeff     = 0.0",
                "couple_moments_momentum_coeff",
            ),
            (
                "crpai_prolate",
                "pic_deltaf_f0                     = kappa_aniso",
                "pic_deltaf_f0                     = kappa_iso",
                "pic_deltaf_f0",
            ),
            (
                "crpai_oblate",
                "blocked_pending_manuscript_text_caption_review",
                "reviewed_right_handed",
                "handedness_mapping",
            ),
        )
        for case, old, new, expected in mutations:
            with self.subTest(case=case, expected=expected):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "candidate.athinput"
                    text = q007.DECKS[case].read_text(encoding="utf-8")
                    self.assertIn(old, text)
                    path.write_text(text.replace(old, new, 1), encoding="utf-8")
                    with self.assertRaisesRegex(q007.ContractError, expected):
                        q007.validate_deck(path, case)

    def test_exact_deck_map_rejects_effectful_optional_controls(self) -> None:
        additions = (
            ("\n<mhd>\nviscosity = 1.0\n", "mhd.*unexpected"),
            ("\n<mhd>\nconst_accel = true\n", "mhd.*unexpected"),
            ("\n<mhd>\nnscalars = 1\n", "mhd.*unexpected"),
            ("\n<shearing_box>\nqshear = 1.5\nomega0 = 1.0\n", "unexpected=.*shearing_box"),
        )
        for addition, expected in additions:
            with self.subTest(addition=addition):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "candidate.athinput"
                    path.write_text(
                        q007.DECKS["crsi"].read_text(encoding="utf-8") + addition,
                        encoding="utf-8",
                    )
                    with self.assertRaisesRegex(q007.ContractError, expected):
                        q007.validate_deck(path, "crsi")

    def test_source_contract_is_additive_and_narrow(self) -> None:
        source = q007.validate_source_contract()
        self.assertTrue(source["compilation_unit_registered"])
        self.assertTrue(source["fresh_and_restart_dispatch_registered"])
        self.assertTrue(source["narrow_exact_isothermal_true_deltaf_parser_allowance"])
        self.assertTrue(source["bounded_serial_replay_guarded"])
        self.assertTrue(source["effectful_optional_mhd_controls_rejected"])

    def test_report_keeps_qualification_boundaries_open(self) -> None:
        report = q007.build_preparation_report()
        boundary = report["nonqualification_boundary"]
        self.assertTrue(boundary["deck_source_freeze"])
        self.assertTrue(boundary["true_deltaf_parser_path"])
        self.assertTrue(boundary["exact_isothermal_mhd_startup"])
        self.assertTrue(boundary["source_local_log_bin_weighted_loading"])
        self.assertTrue(boundary["source_local_deterministic_four_branch_wave_spectrum"])
        self.assertTrue(boundary["bounded_serial_crsi_runtime_replay"])
        self.assertTrue(boundary["full_q1_q2_dispersion_oracle_prepared"])
        for name in (
            "paper_literal_sampling_algorithm_reviewed",
            "paper_literal_wave_seed_and_discrete_mode_set_reviewed",
            "runtime_growth_fit",
            "crpai_handedness_mapping_reviewed",
            "mpi_qualification",
            "gpu_qualification",
            "frontier_authorization",
            "external_review",
            "qualifying_evidence",
        ):
            self.assertFalse(boundary[name])
        self.assertFalse(report["claim_closure"])
        self.assertFalse(report["frontier_authorization"])

    def test_log_bin_weighting_is_normalized_and_deterministic(self) -> None:
        loading = q007.loading_contract(300.0, 1.25)
        self.assertEqual(loading["bin_count"], 8)
        self.assertEqual(loading["particles_per_cell_per_bin"], 256)
        self.assertAlmostEqual(loading["shell_fraction_sum"], 1.0)
        self.assertAlmostEqual(loading["bins"][0]["lower"], 300.0 / 500.0)
        self.assertAlmostEqual(loading["bins"][-1]["upper"], 300.0 * 500.0)
        self.assertTrue(
            all(item["particle_macro_weight"] > 0.0 for item in loading["bins"])
        )
        self.assertFalse(loading["paper_literal_sampling_algorithm_claimed"])

    def test_startup_loading_matches_sparse_independent_pins(self) -> None:
        for case, tag, species, point, state, velocity, weight, f0 in (
            _PINNED_STARTUP_LOADING_SAMPLES
        ):
            with self.subTest(case=case, tag=tag):
                loading = q007._expected_startup_loading(case)
                self.assertEqual(loading["tags"][tag], tag)
                self.assertEqual(loading["species"][tag], species)
                np.testing.assert_allclose(loading["points"][tag], point, rtol=0.0)
                np.testing.assert_allclose(loading["states"][tag], state, rtol=1.0e-13)
                np.testing.assert_allclose(
                    loading["velocities"][tag], velocity, rtol=1.0e-13
                )
                self.assertTrue(math.isclose(loading["weights"][tag], weight,
                                             rel_tol=1.0e-13))
                self.assertTrue(math.isclose(loading["f0"][tag], f0,
                                             rel_tol=1.0e-13))

    def test_four_branch_wave_carriers_recover_discrete_branch_power(self) -> None:
        for length, dx, seed in ((320.0, 10.0, 700702), (640.0, 20.0, 700704)):
            with self.subTest(length=length):
                x1 = (np.arange(32) + 0.5) * dx
                state = np.array([
                    q007.four_branch_alfven_state(x, length, 1.0e-3, seed)
                    for x in x1
                ])
                payload = {
                    "bcc2": state[:, 0],
                    "bcc3": state[:, 1],
                    "vely": state[:, 2],
                    "velz": state[:, 3],
                }
                q007._validate_initial_branch_spectrum(q007._branch_spectrum(payload))
                payload["bcc2"] = payload["bcc2"].copy()
                payload["bcc2"][0] += 1.0e-2
                with self.assertRaisesRegex(q007.ContractError, "spectrum drifted"):
                    q007._validate_initial_branch_spectrum(q007._branch_spectrum(payload))

    def test_initial_mhd_payload_rejects_mutations_hidden_by_branch_power(self) -> None:
        for target in ("bcc2_dc", "bcc2_transverse", "velx"):
            with self.subTest(target=target):
                payload = _startup_mhd_payload("crsi")
                if target == "bcc2_dc":
                    payload["bcc2"] += 1.0e-3
                elif target == "bcc2_transverse":
                    payload["bcc2"][0, 0, 0] += 1.0e-3
                else:
                    payload["velx"][0, 0, 0] += 1.0e-3
                with self.assertRaisesRegex(q007.ContractError, "cells drifted"):
                    q007._validate_initial_mhd_payload(payload, "crsi")

    def test_initial_mhd_payload_rejects_grid_metadata_and_shape_drift(self) -> None:
        mutations = (
            ("MaxLevel", "refinement level"),
            ("x1f", "x1f coordinates"),
            ("x2v", "x2v coordinates"),
            ("shape", "shape drifted"),
        )
        for case in ("crsi", "crpai_prolate", "crpai_oblate"):
            for target, expected in mutations:
                with self.subTest(case=case, target=target):
                    payload = _startup_mhd_payload(case)
                    if target == "MaxLevel":
                        payload["MaxLevel"] = 1
                    elif target == "x1f":
                        payload["x1f"][0] += 1.0
                    elif target == "x2v":
                        payload["x2v"][0] += 1.0
                    else:
                        payload["dens"] = payload["dens"].reshape(4, 1, 32)
                    with self.assertRaisesRegex(q007.ContractError, expected):
                        q007._validate_initial_mhd_payload(payload, case)

    def test_full_q1_q2_oracle_is_finite_and_keeps_labels_open(self) -> None:
        oracle = q007.full_q1_q2_dispersion_oracle()
        self.assertEqual(oracle["mode_count"], 8)
        self.assertFalse(oracle["assigns_crpai_handedness_labels"])
        self.assertFalse(oracle["runtime_growth_fit_claimed"])
        crsi = oracle["tables"]["crsi"]
        crpai = oracle["tables"]["crpai"]
        self.assertEqual(crsi["carrier_length"], 320.0)
        self.assertEqual(crpai["carrier_length"], 640.0)
        self.assertAlmostEqual(crsi["modes"][0]["k"], 2.0 * math.pi / 320.0)
        self.assertAlmostEqual(crpai["modes"][0]["k"], 2.0 * math.pi / 640.0)
        for mode in crsi["modes"]:
            self.assertGreater(mode["q1_crsi"], 0.0)
            self.assertGreater(mode["q2_crsi"], 0.0)
            for polarization in ("-1", "1"):
                for direction in ("forward", "backward"):
                    root = mode["crsi_signed_polarizations"][polarization][direction]
                    self.assertTrue(math.isfinite(root["real"]))
                    self.assertTrue(math.isfinite(root["imag"]))
        for mode in crpai["modes"]:
            self.assertGreater(mode["q1_crpai"], 0.0)
            self.assertGreater(mode["q2_crpai"], 0.0)
            for role in ("prolate", "oblate"):
                for polarization in ("-1", "1"):
                    for direction in ("forward", "backward"):
                        root = mode["crpai_signed_polarizations"][role][polarization][
                            direction
                        ]
                        self.assertTrue(math.isfinite(root["real"]))
                        self.assertTrue(math.isfinite(root["imag"]))

    def test_startup_particle_payload_validates_all_three_cases(self) -> None:
        for case in ("crsi", "crpai_prolate", "crpai_oblate"):
            with self.subTest(case=case):
                report = q007._validate_startup_particle_payload(
                    _startup_particle_payload(case), case
                )
                self.assertEqual(report["point_layout_max_abs_error"], 0.0)
                self.assertEqual(report["deltaf_f0_max_abs_error"], 0.0)
                self.assertEqual(report["vtk_serialized_velocity_max_abs_error"], 0.0)
                self.assertEqual(
                    report["transformed_shell_serialization_max_abs_error"], 0.0
                )
                self.assertLess(report["antipodal_velocity_pair_sum_max_abs"], 1.0e-8)

    def test_startup_particle_payload_mutations_fail_closed(self) -> None:
        mutations = (
            ("points", "center-distribution points drifted"),
            ("vel", "deterministic angular sampler drifted"),
            ("deltaf_f0", "delta-f f0 shape drifted"),
        )
        for target, expected in mutations:
            with self.subTest(target=target):
                payload = _startup_particle_payload("crsi")
                if target == "points":
                    payload.points[0, 0] += 1.0
                elif target == "vel":
                    payload.vectors["vel"][0, 0] += 1.0
                else:
                    payload.scalars["deltaf_f0"][0] *= 0.5
                with self.assertRaisesRegex(q007.ContractError, expected):
                    q007._validate_startup_particle_payload(payload, "crsi")

    def test_particle_payload_rejects_round_robin_species_drift(self) -> None:
        payload = _startup_particle_payload("crsi")
        payload.scalars["species"][0] = 1
        with patch.object(q007, "_read_particle_vtk", return_value=payload):
            with self.assertRaisesRegex(q007.ContractError, "round-robin"):
                q007._particle_payload(
                    Path("unused"),
                    case="crsi",
                    require_initial_df_zero=True,
                    validate_startup_loading=True,
                )

    def test_runtime_extractor_rejects_non_orion_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(q007.ContractError, "remain below"):
                q007._authorized_orion_runtime_root(Path(directory))

    def test_runtime_tree_freeze_and_verify_fail_closed(self) -> None:
        original_root = q007.ORION_PIC_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "runtime"
            tree.mkdir()
            (tree / "payload").mkdir()
            (tree / "payload/result.json").write_text("{}\n", encoding="utf-8")
            q007.ORION_PIC_ROOT = Path(directory)
            try:
                frozen = q007.freeze_runtime_tree(tree)
                self.assertTrue(frozen["recursively_read_only"])
                self.assertEqual(frozen["inventoried_file_count"], 2)
                self.assertEqual(frozen["writable_entries"], [])
                verified = q007.verify_frozen_runtime_tree(
                    tree, frozen["inventory_sha256"]
                )
                self.assertEqual(verified["inventory_sha256"], frozen["inventory_sha256"])

                inventory = tree / q007.INVENTORY_NAME
                inventory.chmod(inventory.stat().st_mode | stat.S_IWUSR)
                inventory.write_text(
                    inventory.read_text(encoding="utf-8")
                    + "0" * 64
                    + "  nested/"
                    + q007.INVENTORY_NAME
                    + "\n",
                    encoding="utf-8",
                )
                inventory.chmod(inventory.stat().st_mode & ~stat.S_IWUSR)
                with self.assertRaisesRegex(q007.ContractError, "not sorted|membership drifted"):
                    q007.verify_frozen_runtime_tree(tree, _sha256(inventory))
            finally:
                for path in [tree, *tree.rglob("*")]:
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
                q007.ORION_PIC_ROOT = original_root

    def test_runtime_tree_verify_rejects_symlink_and_unsafe_inventory_path(self) -> None:
        original_root = q007.ORION_PIC_ROOT
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory) / "runtime"
            tree.mkdir()
            (tree / "payload").write_text("bounded\n", encoding="utf-8")
            q007.ORION_PIC_ROOT = Path(directory)
            try:
                frozen = q007.freeze_runtime_tree(tree)
                tree.chmod(tree.stat().st_mode | stat.S_IWUSR)
                inventory = tree / q007.INVENTORY_NAME
                inventory.chmod(inventory.stat().st_mode | stat.S_IWUSR)
                original_inventory = inventory.read_text(encoding="utf-8")
                inventory.write_text(
                    "0" * 64 + "  ../outside\n", encoding="utf-8"
                )
                inventory.chmod(inventory.stat().st_mode & ~stat.S_IWUSR)
                tree.chmod(tree.stat().st_mode & ~stat.S_IWUSR)
                with self.assertRaisesRegex(q007.ContractError, "unsafe"):
                    q007.verify_frozen_runtime_tree(tree, _sha256(inventory))

                tree.chmod(tree.stat().st_mode | stat.S_IWUSR)
                inventory.chmod(inventory.stat().st_mode | stat.S_IWUSR)
                inventory.write_text(original_inventory, encoding="utf-8")
                inventory.chmod(inventory.stat().st_mode & ~stat.S_IWUSR)
                (tree / "alias").symlink_to(tree / "payload")
                tree.chmod(tree.stat().st_mode & ~stat.S_IWUSR)
                with self.assertRaisesRegex(q007.ContractError, "symlink"):
                    q007.verify_frozen_runtime_tree(tree, frozen["inventory_sha256"])
            finally:
                for path in [tree, *tree.rglob("*")]:
                    if not path.is_symlink():
                        path.chmod(path.stat().st_mode | stat.S_IWUSR)
                q007.ORION_PIC_ROOT = original_root

    def test_runtime_invocation_closure_rejects_manifest_mutations(self) -> None:
        mutations = (
            ("crsi", "executable_sha256", "3" * 64, "executable_sha256"),
            ("crsi", "overrides", ["time/nlim=1"], "overrides"),
            ("crsi", "returncode", 1, "returncode"),
            ("crsi", "returncode", False, "returncode"),
            ("crsi", "returncode", 0.0, "returncode"),
            ("crsi", "produced_file_sha256", {}, "produced payload"),
        )
        for label, key, value, expected in mutations:
            with self.subTest(label=label, key=key):
                with tempfile.TemporaryDirectory() as directory:
                    root = Path(directory) / "runtime"
                    pinned = _write_runtime_invocation_tree(root)
                    q007._validate_runtime_invocations(root, pinned)
                    invocation_path = root / label / "invocation.json"
                    invocation = json.loads(invocation_path.read_text(encoding="utf-8"))
                    invocation[key] = value
                    invocation_path.write_text(
                        json.dumps(invocation, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8",
                    )
                    with self.assertRaisesRegex(q007.ContractError, expected):
                        q007._validate_runtime_invocations(root, pinned)

    def test_retained_json_guard_rejects_boolean_aliases(self) -> None:
        expected = {
            "clean_empty_directory_configure": True,
            "configure_returncode": 0,
            "mpi_used": False,
        }
        for replacement in (
            {"clean_empty_directory_configure": 1},
            {"clean_empty_directory_configure": 1.0},
            {"configure_returncode": False},
            {"configure_returncode": 0.0},
            {"mpi_used": 0},
            {"mpi_used": 0.0},
        ):
            with self.subTest(replacement=replacement):
                with self.assertRaisesRegex(q007.ContractError, "primitive type drifted"):
                    q007._require_exact_match({**expected, **replacement}, expected, "metadata")

    def test_runtime_invocation_closure_rejects_unexpected_claim_key(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "runtime"
            pinned = _write_runtime_invocation_tree(root)
            invocation_path = root / "crsi/invocation.json"
            invocation = json.loads(invocation_path.read_text(encoding="utf-8"))
            invocation["qualifying_evidence"] = True
            invocation_path.write_text(
                json.dumps(invocation, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(q007.ContractError, "schema drifted"):
                q007._validate_runtime_invocations(root, pinned)

    def test_runtime_invocation_closure_rejects_retained_payload_mutations(self) -> None:
        for target, expected in (
            ("returncode", "returncode text"),
            ("unexpected", "payload topology"),
            ("deck", "retained deck"),
        ):
            with self.subTest(target=target):
                with tempfile.TemporaryDirectory() as directory:
                    root = Path(directory) / "runtime"
                    pinned = _write_runtime_invocation_tree(root)
                    if target == "returncode":
                        (root / "crsi/returncode.txt").write_text("1\n", encoding="utf-8")
                    elif target == "unexpected":
                        (root / "crsi/unexpected.txt").write_text("extra\n", encoding="utf-8")
                    else:
                        deck = root / "decks/pic_q007_paper_crsi_linear_preparation.athinput"
                        deck.write_text(
                            deck.read_text(encoding="utf-8").replace(
                                "nlim       = 2", "nlim       = 3", 1
                            ),
                            encoding="utf-8",
                        )
                        invocation = root / "crsi/invocation.json"
                        payload = json.loads(invocation.read_text(encoding="utf-8"))
                        payload["deck"]["sha256"] = _sha256(deck)
                        invocation.write_text(
                            json.dumps(payload, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8",
                        )
                    with self.assertRaisesRegex(q007.ContractError, expected):
                        q007._validate_runtime_invocations(root, pinned)

    def test_runtime_tree_topology_rejects_unexpected_empty_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "runtime"
            _write_runtime_invocation_tree(root)
            for name in (
                q007.INVENTORY_NAME,
                q007.FREEZE_RECEIPT_NAME,
                q007.PINNED_EXECUTABLE_BINDING_NAME,
            ):
                (root / name).write_text("bounded\n", encoding="utf-8")
            q007._validate_runtime_tree_topology(root)
            (root / "unexpected-empty").mkdir()
            with self.assertRaisesRegex(q007.ContractError, "directory topology"):
                q007._validate_runtime_tree_topology(root)

    def test_runtime_extraction_cli_requires_anchored_inventory_digest(self) -> None:
        with patch.object(sys, "argv", [str(q007.__file__), "--runtime-root", "/tmp/runtime"]):
            with self.assertRaisesRegex(q007.ContractError, "requires an anchored digest"):
                q007.main()

    def test_strict_parser_rejects_duplicate_parameter(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.athinput"
            path.write_text("<mesh>\nnx1 = 32\nnx1 = 32\n", encoding="utf-8")
            with self.assertRaisesRegex(q007.ContractError, "duplicate parameter"):
                q007.parse_athinput(path)

    def test_readiness_sidecar_binds_only_new_q007_artifacts(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-007")
        self.assertEqual(sidecar["qualification_effect"], q007.QUALIFICATION_EFFECT)
        self.assertFalse(sidecar["claim_closure"])
        self.assertEqual(sidecar["frontier_authorization"], "not_bound")
        expected_paths = {
            "src/pgen/tests/q007_paper_deltaf_linear.hpp",
            "src/pgen/tests/q007_paper_deltaf_linear.cpp",
            "inputs/tests/pic_q007_paper_crsi_linear_preparation.athinput",
            "inputs/tests/pic_q007_paper_crpai_linear_prolate_preparation.athinput",
            "inputs/tests/pic_q007_paper_crpai_linear_oblate_preparation.athinput",
            "tst/publication/immutable_orion_tree.py",
            "tst/publication/analyze_q007_paper_deltaf_linear_preparation.py",
            "tst/publication/test_analyze_q007_paper_deltaf_linear_preparation.py",
        }
        bindings = sidecar["artifact_bindings"]
        self.assertEqual(set(bindings), expected_paths)
        for relative, expected in bindings.items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected)
        shared_paths = {
            "docs/source/engineering/pic_mhd_model_contract.md",
            "src/CMakeLists.txt",
            "src/particles/particles.cpp",
            "src/pgen/pgen.cpp",
            "src/pgen/pgen.hpp",
            "src/outputs/vtk_prtcl.cpp",
            "tst/publication/pvtk_particles.py",
        }
        self.assertEqual(set(sidecar["shared_artifact_bindings"]), shared_paths)
        for relative, expected in sidecar["shared_artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected)
        boundary = sidecar["nonqualification_boundary"]
        self.assertFalse(boundary["qualifying_evidence"])
        self.assertTrue(boundary["bounded_serial_crsi_runtime_replay"])
        self.assertTrue(boundary["source_local_log_bin_weighted_loading"])
        self.assertTrue(boundary["source_local_deterministic_four_branch_wave_spectrum"])
        self.assertFalse(boundary["runtime_growth_fit"])


if __name__ == "__main__":
    unittest.main()
