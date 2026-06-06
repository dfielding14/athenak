"""Focused tests for deterministic Stage I scientific products."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_scientific_products.py"
ACCEPTANCE_UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_scientific_acceptance.py"
PAPER_ANALYZER = REPOSITORY / "scripts/analyze_cgl_lf_paper.py"
BINARY_PARSER = REPOSITORY / "vis/python/bin_convert.py"


def load_utility(path: Path = UTILITY, name: str = "cgl_lf_stage_i_scientific_products"):
    """Load one standalone Stage I utility without importing the dirty analyzer."""

    spec = importlib.util.spec_from_file_location(
        name, path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


products = load_utility()
acceptance = load_utility(
    ACCEPTANCE_UTILITY, "cgl_lf_stage_i_scientific_acceptance_for_products"
)
paper_analyzer = load_utility(PAPER_ANALYZER, "analyze_cgl_lf_paper_for_products")


def sha256(path: Path) -> str:
    """Return one file digest."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    """Write deterministic JSON fixture bytes."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_history(path: Path, columns: dict[str, list[float]]) -> None:
    """Write one minimal Athena history fixture."""

    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    lines = [
        "# " + " ".join(
            f"[{index}]={name}" for index, name in enumerate(names, start=1)
        )
    ]
    lines.extend(
        " ".join(format(value, ".17g") for value in row)
        for row in zip(*(columns[name] for name in names))
    )
    path.write_text("\n".join(lines) + "\n")


def fixture_histories(
    root: Path,
) -> tuple[Path, Path, dict[str, list[float]], dict[str, list[float]]]:
    """Create histories with all required scientific-product columns."""

    times = [0.0, 0.25, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
    count = len(times)
    user = {
        "time": times,
        "volume": [2.0] * count,
        "abs_dp": [0.2 + 0.01 * index for index in range(count)],
        "beta": [20.0] * count,
        "fire_vol": [0.004] * count,
        "force_pwr": [0.64] * count,
        "hard_vol": [0.0] * count,
        "kinetic": [2.0] * count,
        "magnetic": [3.0] * count,
        "mirror_vol": [0.002] * count,
        "nu_eff": [20.0] * count,
        "force_prp2": [1.0] * count,
        "force_prl2": [0.0] * count,
    }
    mhd = {
        "time": times,
        "lf_hwproj": [float(index) for index in range(count)],
        "lf_cpwrk": [1.0e-4 * index for index in range(count)],
        "lf_cawrk": [-2.0e-4 * index for index in range(count)],
        "lf_qface": [100.0 * index for index in range(count)],
        "lf_qprwrk": [1.0e-3 * index for index in range(count)],
        "lf_qpewrk": [-2.0e-3 * index for index in range(count)],
    }
    mhd_path = root / "bundle/history/case.mhd.hst"
    user_path = root / "bundle/history/case.user.hst"
    write_history(mhd_path, mhd)
    write_history(user_path, user)
    return mhd_path, user_path, mhd, user


def write_curve(path: Path, rows: list[tuple[float, float, float]]) -> None:
    """Write one x/y/uncertainty reference curve."""

    path.parent.mkdir(parents=True, exist_ok=True)
    text = "x,y,y_uncertainty\n"
    text += "\n".join(",".join(format(value, ".17g") for value in row) for row in rows)
    path.write_text(text + "\n")


def build_fixture(root: Path) -> tuple[dict[str, object], dict[str, Path]]:
    """Build an internally authenticated offline accepted bundle."""

    mhd_path, user_path, mhd, user = fixture_histories(root)
    references = root / "references"
    unstable_csv = references / "unstable.csv"
    density_csv = references / "density.csv"
    write_curve(unstable_csv, [(8.0, 0.003, 0.001), (9.0, 0.003, 0.001), (10.0, 0.003, 0.001)])
    write_curve(density_csv, [(-0.1, 1.0, 0.5), (0.0, 2.0, 0.5), (0.1, 1.0, 0.5)])
    reference_manifest = references / "curves.json"
    write_json(reference_manifest, {
        "schema_version": 1,
        "provenance": {
            "method": "digitized",
            "source_description": "fixture",
            "uncertainty_description": "fixture",
            "digitization_tool": "fixture",
            "source_figure": "fixture.pdf",
            "source_figure_sha256": "a" * 64,
        },
        "curves": [
            {
                "id": "fixture_unstable",
                "case": "fixture_case",
                "product": "history.unstable_fraction",
                "data_file": unstable_csv.name,
                "data_sha256": sha256(unstable_csv),
                "interpolation": "linear",
            },
            {
                "id": "fixture_density",
                "case": "fixture_case",
                "product": "pdf.density_fluctuation",
                "data_file": density_csv.name,
                "data_sha256": sha256(density_csv),
                "interpolation": "linear",
            },
        ],
        "surfaces": [],
    })
    stage_manifest = root / "stage_manifest.json"
    write_json(stage_manifest, {
        "cases": [
            {
                "id": "R02",
                "name": "fixture_case",
                "input": "inputs/fixture.athinput",
            }
        ],
        "panel_status": {
            "analysis_case_aliases": {},
            "reference_manifests": {
                "fixture": {
                    "path": reference_manifest.name,
                    "sha256": sha256(reference_manifest),
                }
            },
            "reference_product_bindings": {
                "fixture_unstable": {
                    "kind": "curve",
                    "case": "fixture_case",
                    "product": "history.unstable_fraction",
                    "data_file": unstable_csv.name,
                    "data_sha256": sha256(unstable_csv),
                    "reference_manifest": "fixture",
                },
                "fixture_density": {
                    "kind": "curve",
                    "case": "fixture_case",
                    "product": "pdf.density_fluctuation",
                    "data_file": density_csv.name,
                    "data_sha256": sha256(density_csv),
                    "reference_manifest": "fixture",
                },
            },
            "panels": [
                {"id": "fig2b", "reference_products": ["fixture_unstable"]},
                {"id": "fig4a", "reference_products": ["fixture_density"]},
            ],
        },
    })
    stage_sha = sha256(stage_manifest)
    segments = [root / "segments/s00.json", root / "segments/s01.json"]
    segment_indices = ([0, 1], list(range(1, len(mhd["time"]))))
    segment_histories: list[tuple[Path, Path]] = []
    for index, indices in enumerate(segment_indices):
        mhd_segment = root / f"segment_output/s0{index}.mhd.hst"
        user_segment = root / f"segment_output/s0{index}.user.hst"
        write_history(
            mhd_segment,
            {name: [values[row] for row in indices] for name, values in mhd.items()},
        )
        write_history(
            user_segment,
            {name: [values[row] for row in indices] for name, values in user.items()},
        )
        segment_histories.append((mhd_segment, user_segment))
    executable_sha = "b" * 64
    input_sha = "c" * 64
    for index, (path, final_time) in enumerate(zip(segments, (0.25, 10.0))):
        parent = None
        if index:
            parent = {
                "manifest": str(segments[index - 1].resolve()),
                "case_id": "R02",
                "result": "accepted",
                "final_time": 0.25,
                "executable_sha256": executable_sha,
                "input_sha256": input_sha,
            }
        write_json(path, {
            "accounting": {
                "result": "accepted",
                "case_id": "R02",
                "case_name": "fixture_case",
                "segment": f"s0{index}_fixture",
                "executable_sha256": executable_sha,
            },
            "command": {
                "executable_sha256": executable_sha,
                "input_sha256": input_sha,
                "matrix_sha256": stage_sha,
                "parent_segment": parent,
                "restart_file": None if index == 0 else "parent.rst",
                "restart_files": [] if index == 0 else ["parent.rst"],
            },
            "scientific_inspection": {
                "accepted": True,
                "case_id": "R02",
                "manifest": str(path.resolve()),
                "final_time": final_time,
                "mhd_history": {
                    "path": str(segment_histories[index][0].resolve()),
                    "sha256": sha256(segment_histories[index][0]),
                    "size_bytes": segment_histories[index][0].stat().st_size,
                },
                "user_history": {
                    "path": str(segment_histories[index][1].resolve()),
                    "sha256": sha256(segment_histories[index][1]),
                    "size_bytes": segment_histories[index][1].stat().st_size,
                },
                "snapshot_times": [],
                "snapshots": [],
            },
        })
    bundle = root / "bundle/manifest.json"
    write_json(bundle, {
        "workflow": "paper-mks24-stage-i-production",
        "status": "accepted_for_analysis",
        "production_case_id": "R02",
        "required_final_time": 10.0,
        "accepted_final_time": 10.0,
        "stage_i_manifest": str(stage_manifest.resolve()),
        "cases": [
            {
                "name": "fixture_case",
                "input": "inputs/fixture.athinput",
                "status": "passed",
                "model_choices": {"forcing_tcorr": "2.0"},
                "outputs": {
                    "mhd_history": "history/case.mhd.hst",
                    "user_history": "history/case.user.hst",
                    "snapshot_paths": [],
                },
            }
        ],
        "production_segment_manifests": [str(path.resolve()) for path in segments],
    })
    request = {
        "authority_mode": "offline",
        "bundle_manifest": str(bundle.resolve()),
        "expected_bundle_sha256": sha256(bundle),
        "time_start": 8.0,
        "time_end": 10.0,
        "reference_root": str(references.resolve()),
        "snapshot_mode": "skip",
        "pdf_bins": 16,
        "alignment_shells": [1, 2, 4],
        "max_snapshot_cells": 1_000_000,
    }
    return request, {
        "bundle": bundle,
        "mhd": mhd_path,
        "user": user_path,
        "segment_user": segment_histories[1][1],
        "unstable_csv": unstable_csv,
        "density_csv": density_csv,
    }


def retained_evidence(tmp_path: Path) -> tuple[dict[str, object], Path, dict[str, Path]]:
    """Generate and retain one fixture evidence record."""

    request, paths = build_fixture(tmp_path)
    evidence = products.build_evidence(request)
    path = tmp_path / "evidence.json"
    products.write_candidate(path, evidence)
    return evidence, path, paths


def test_generator_does_not_import_dirty_analyzer():
    """The standalone utility must not import or mention analyzer internals."""

    source = UTILITY.read_text()
    assert "import analyze_cgl_lf_paper" not in source
    assert "from analyze_cgl_lf_paper" not in source


def test_athena_binary_parser_is_authority_bound_before_execution(tmp_path):
    """A replacement parser is rejected before any of its bytes execute."""

    copied_utility = (
        tmp_path / "repo/scripts/frontier/cgl_lf_stage_i_scientific_products.py"
    )
    copied_parser = tmp_path / "repo/vis/python/bin_convert.py"
    copied_utility.parent.mkdir(parents=True)
    copied_parser.parent.mkdir(parents=True)
    copied_utility.write_bytes(UTILITY.read_bytes())
    marker = tmp_path / "untrusted-parser-executed"
    copied_parser.write_text(
        "from pathlib import Path\n"
        f"Path({str(marker)!r}).write_text('executed')\n"
    )

    with pytest.raises(RuntimeError, match="Athena binary parser SHA-256 differs"):
        load_utility(copied_utility, "cgl_lf_stage_i_scientific_products_replaced_parser")
    assert marker.exists() is False


def test_generated_evidence_binds_exact_athena_binary_parser(tmp_path):
    """Every generated record binds the exact parser loaded from reviewed bytes."""

    evidence, _, _ = retained_evidence(tmp_path)
    assert evidence["runtime"]["binary_parser"] == products.BIN_CONVERT_BINDING
    assert evidence["runtime"]["binary_parser"]["sha256"] == sha256(BINARY_PARSER)


def test_history_products_reference_residuals_and_deferred_products(tmp_path):
    """History products are recomputed while unsupported snapshot products fail closed."""

    evidence, _, _ = retained_evidence(tmp_path)
    metrics = evidence["scientific_acceptance_metrics"]
    assert metrics["kinetic"]["time_weighted_mean"] == pytest.approx(2.0)
    assert metrics["unstable_occupancy"]["time_weighted_mean"] == pytest.approx(0.003)
    comparisons = evidence["reference_curve_comparisons"]["comparisons"]
    assert comparisons["fixture_unstable"]["available"] is True
    assert comparisons["fixture_unstable"]["normalized_residual_rms"] == pytest.approx(0.0)
    assert comparisons["fixture_density"]["available"] is False
    assert evidence["snapshot_ensemble"]["result"] == "inconclusive"
    assert evidence["result"] == "inconclusive"


def test_generated_evidence_emits_nested_schema_two_acceptance_contract(tmp_path):
    """Generated products directly expose the contract and case shape acceptance consumes."""

    evidence, _, paths = retained_evidence(tmp_path)
    contract = evidence["scientific_acceptance_contract"]
    assert contract == {
        "schema_version": 2,
        "record_type": "stage-i-scientific-acceptance-reviewed-products",
        "case_id": "R02",
        "case_name": "fixture_case",
        "analysis_window": {"time_start": 8.0, "time_end": 10.0},
        "accepted_bundle_manifest": evidence["authentication"]["bundle_manifest"],
        "generator": evidence["generator"],
        "deterministic_replay_verification": (
            "exact-semantic-replay-from-bound-canonical-case-inputs"
        ),
        "stage_i_manifest_sha256": evidence["authentication"]["stage_i_manifest"]["sha256"],
        "mhd_history": {
            "path": str(paths["mhd"].resolve()),
            "size_bytes": paths["mhd"].stat().st_size,
            "sha256": sha256(paths["mhd"]),
        },
        "user_history": {
            "path": str(paths["user"].resolve()),
            "size_bytes": paths["user"].stat().st_size,
            "sha256": sha256(paths["user"]),
        },
    }
    selected = evidence["cases"]["fixture_case"]
    assert selected["analysis_window"] == contract["analysis_window"]
    assert selected["lf_counter_increments"] == evidence["lf_counter_increments"]
    assert selected["snapshot_ensemble"] == evidence["snapshot_ensemble"]
    assert (
        selected["scientific_acceptance_convergence"]
        == evidence["scientific_acceptance_convergence"]
    )


def test_replay_is_byte_identical_and_deterministic(tmp_path):
    """Replay recomputes the record and compares exact serialized bytes."""

    evidence, path, _ = retained_evidence(tmp_path)
    expected = sha256(path)
    replayed, identical = products.replay_evidence(path, expected)
    assert identical is True
    assert replayed == evidence
    assert products.stable_json(products.build_evidence(evidence["request"])) == path.read_bytes()


def test_forged_nested_acceptance_contract_cannot_pass_replay(tmp_path):
    """A self-digested forged schema-2 contract still fails deterministic replay."""

    evidence, _, _ = retained_evidence(tmp_path)
    forged = deepcopy(evidence)
    forged.pop("evidence_digest")
    forged["scientific_acceptance_contract"]["analysis_window"]["time_start"] = 7.0
    forged = products.seal_evidence(forged)
    path = tmp_path / "forged-contract.json"
    products.write_candidate(path, forged)
    _, identical = products.replay_evidence(path, sha256(path))
    assert identical is False


def test_forged_claimed_output_cannot_pass_replay(tmp_path):
    """A self-digested hand-authored result still fails full recomputation."""

    evidence, _, _ = retained_evidence(tmp_path)
    forged = deepcopy(evidence)
    forged.pop("evidence_digest")
    forged["scientific_acceptance_metrics"]["kinetic"]["time_weighted_mean"] = 999.0
    forged = products.seal_evidence(forged)
    path = tmp_path / "forged.json"
    products.write_candidate(path, forged)
    _, identical = products.replay_evidence(path, sha256(path))
    assert identical is False


def test_forged_request_cannot_pass_replay(tmp_path):
    """A forged request is rejected even when its evidence self-digest is repaired."""

    evidence, _, _ = retained_evidence(tmp_path)
    forged = deepcopy(evidence)
    forged.pop("evidence_digest")
    forged["request"]["expected_bundle_sha256"] = "0" * 64
    forged = products.seal_evidence(forged)
    path = tmp_path / "forged-request.json"
    products.write_candidate(path, forged)
    with pytest.raises(products.ScientificProductsError, match="whole-case bundle manifest SHA"):
        products.replay_evidence(path, sha256(path))


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_replay_request_rejects_missing_or_extra_fields(tmp_path, mutation):
    """Replay requests are closed to the exact reviewed key set."""

    evidence, _, _ = retained_evidence(tmp_path)
    forged = deepcopy(evidence)
    forged.pop("evidence_digest")
    if mutation == "missing":
        forged["request"].pop("snapshot_mode")
    else:
        forged["request"]["unreviewed_extension"] = True
    path = tmp_path / f"forged-{mutation}-request.json"
    products.write_candidate(path, products.seal_evidence(forged))
    with pytest.raises(products.ScientificProductsError, match="request keys differ"):
        products.replay_evidence(path, sha256(path))


def test_bundle_history_input_mutation_is_rejected(tmp_path):
    """Mutating a merged bundle history breaks its authenticated segment merge."""

    _, path, paths = retained_evidence(tmp_path)
    expected = sha256(path)
    paths["user"].write_text(paths["user"].read_text().replace("20 ", "21 ", 1))
    with pytest.raises(products.ScientificProductsError, match="user history differs"):
        products.replay_evidence(path, expected)


def test_segment_history_input_mutation_is_rejected(tmp_path):
    """Mutating one accepted segment history breaks its inspection digest."""

    _, path, paths = retained_evidence(tmp_path)
    expected = sha256(path)
    paths["segment_user"].write_text(
        paths["segment_user"].read_text().replace("20 ", "21 ", 1)
    )
    with pytest.raises(products.ScientificProductsError, match="segment 1 user_history SHA"):
        products.replay_evidence(path, expected)


def test_reference_csv_mutation_is_rejected(tmp_path):
    """Reference CSV bytes are authenticated against the Stage I binding."""

    _, path, paths = retained_evidence(tmp_path)
    expected = sha256(path)
    paths["unstable_csv"].write_text(
        paths["unstable_csv"].read_text().replace("0.003", "0.004", 1)
    )
    with pytest.raises(products.ScientificProductsError, match="reference CSV fixture_unstable SHA"):
        products.replay_evidence(path, expected)


def synthetic_fields(shape: tuple[int, int, int] = (16, 16, 16)) -> dict[str, np.ndarray]:
    """Return one finite periodic CGL snapshot field set."""

    z, y, x = np.meshgrid(
        (np.arange(shape[0]) + 0.5) * 2.0 / shape[0],
        (np.arange(shape[1]) + 0.5) / shape[1],
        (np.arange(shape[2]) + 0.5) / shape[2],
        indexing="ij",
    )
    return {
        "dens": 1.0 + 0.05 * np.sin(2.0 * math.pi * x),
        "velx": np.sin(4.0 * math.pi * x),
        "vely": np.cos(4.0 * math.pi * y),
        "velz": np.sin(2.0 * math.pi * z),
        "eint": 10.0 + 0.1 * np.cos(2.0 * math.pi * z),
        "p_perp": 10.0 + 0.1 * np.sin(2.0 * math.pi * z),
        "bcc1": 0.1 * np.sin(2.0 * math.pi * y),
        "bcc2": 0.1 * np.cos(2.0 * math.pi * x),
        "bcc3": np.ones(shape),
    }


def synthetic_eddy_fields(
    shape: tuple[int, int, int] = (16, 16, 16),
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Return the analyzer's finite local-field eddy synthetic vector fields."""

    z, y, x = np.meshgrid(
        (np.arange(shape[0]) + 0.5) * 2.0 / shape[0],
        (np.arange(shape[1]) + 0.5) / shape[1],
        (np.arange(shape[2]) + 0.5) / shape[2],
        indexing="ij",
    )
    velocity = [
        np.sin(2.0 * math.pi * x) + np.sin(math.pi * z),
        np.sin(2.0 * math.pi * y) + np.sin(math.pi * z),
        np.zeros(shape),
    ]
    magnetic = [
        0.1 * velocity[0],
        0.1 * velocity[1],
        np.ones(shape),
    ]
    return velocity, magnetic


def test_pressure_transfer_has_exact_paper_analyzer_semantic_parity():
    """Stage I preserves the paper analyzer's transfer values and reference curve."""

    fields = synthetic_fields()
    density = fields["dens"]
    velocity = [fields["velx"], fields["vely"], fields["velz"]]
    magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
    delta_p = fields["p_perp"] - fields["eint"]
    lengths = (1.0, 1.0, 2.0)
    dk = math.pi
    expected = paper_analyzer.pressure_transfer(
        density, velocity, magnetic, delta_p, lengths, dk
    )
    actual = products.pressure_transfer(
        density, velocity, magnetic, delta_p, lengths, dk
    )
    for key in (
        "dk",
        "k_perp",
        "transfer",
        "normalization_available",
        "normalization_definition",
        "kinetic_energy",
        "velocity_rms",
        "perpendicular_outer_scale",
        "total_transfer_rate",
        "transfer_normalized_by_total",
        "direct_real_space",
        "shell_sum",
        "closure_error",
    ):
        if isinstance(expected[key], list):
            if key == "k_perp":
                np.testing.assert_array_equal(actual[key], expected[key])
            else:
                np.testing.assert_allclose(
                    actual[key], expected[key], rtol=1.0e-12, atol=2.0e-15
                )
        elif key in ("shell_sum", "closure_error"):
            assert actual[key] == pytest.approx(expected[key], rel=1.0e-12, abs=2.0e-15)
        else:
            assert actual[key] == expected[key]

    actual_ensemble = {"result": "pass", "pressure_transfer": actual}
    expected_ensemble = {"pressure_transfer": expected}
    context = products.BundleContext.__new__(products.BundleContext)
    for product in (
        "pressure_transfer.transfer",
        "pressure_transfer.transfer_normalized_by_total",
    ):
        expected_x, expected_y = paper_analyzer.analyzed_product_curve(
            expected_ensemble, product
        )
        actual_x, actual_y = products.curve_from_product(
            context, actual_ensemble, product
        )
        np.testing.assert_array_equal(actual_x, expected_x)
        np.testing.assert_allclose(actual_y, expected_y, rtol=1.0e-12, atol=2.0e-15)
        reference_x = np.linspace(expected_x[0], expected_x[-1], 11)
        np.testing.assert_allclose(
            products.interpolate_curve(reference_x, actual_x, actual_y, "linear"),
            paper_analyzer.interpolate_analysis_curve(
                reference_x, expected_x, expected_y, "linear"
            ),
            rtol=1.0e-12,
            atol=2.0e-15,
        )


def test_local_field_eddy_anisotropy_has_exact_paper_analyzer_semantic_parity():
    """Stage I preserves paper-analyzer eddy sampling, scales, and reference curves."""

    velocity, magnetic = synthetic_eddy_fields()
    lengths = (1.0, 1.0, 2.0)
    expected = paper_analyzer.local_field_eddy_anisotropy(
        velocity, magnetic, lengths, 200_000, 10, 731
    )
    actual = products.local_field_eddy_anisotropy(
        velocity, magnetic, lengths, 200_000, 10, 731
    )
    for key in (
        "computed",
        "available",
        "definition",
        "conditioning",
        "separation_coordinate",
        "normalization_definition",
        "lperp",
        "angle_degrees",
        "samples_requested",
        "samples_retained",
        "seed",
        "bins",
        "bin_centers_over_lperp",
    ):
        if isinstance(expected[key], list):
            np.testing.assert_array_equal(actual[key], expected[key])
        else:
            assert actual[key] == expected[key]
    for field in ("velocity_perp", "magnetic_perp"):
        for key in (
            "available",
            "ell_perp_over_lperp",
            "ell_parallel_over_lperp",
            "perpendicular_structure_function",
            "parallel_structure_function",
            "perpendicular_sample_counts",
            "parallel_sample_counts",
        ):
            if isinstance(expected[field][key], list):
                np.testing.assert_array_equal(actual[field][key], expected[field][key])
            else:
                assert actual[field][key] == expected[field][key]

    actual_ensemble = {"result": "pass", "eddy_anisotropy": actual}
    expected_ensemble = {"eddy_anisotropy": expected}
    context = products.BundleContext.__new__(products.BundleContext)
    for field in ("velocity_perp", "magnetic_perp"):
        product = f"eddy_anisotropy.{field}"
        expected_x, expected_y = paper_analyzer.analyzed_product_curve(
            expected_ensemble, product
        )
        actual_x, actual_y = products.curve_from_product(
            context, actual_ensemble, product
        )
        np.testing.assert_array_equal(actual_x, expected_x)
        np.testing.assert_array_equal(actual_y, expected_y)
        reference_x = np.linspace(expected_x[0], expected_x[-1], 11)
        np.testing.assert_array_equal(
            products.interpolate_curve(reference_x, actual_x, actual_y, "linear"),
            paper_analyzer.interpolate_analysis_curve(
                reference_x, expected_x, expected_y, "linear"
            ),
        )


def test_pressure_transfer_known_sign_shell_closure_and_product_curve():
    """The MKS24 CGL stress transfer has known sign and closes over k_perp."""

    shape = (16, 16, 16)
    z = np.broadcast_to(
        ((np.arange(shape[0]) + 0.5) * 2.0 / shape[0])[:, None, None],
        shape,
    )
    zero = np.zeros(shape)
    density = np.ones(shape)
    velocity = [zero, zero, np.cos(math.pi * z)]
    magnetic = [zero, zero, np.ones(shape)]
    transfer = products.pressure_transfer(
        density,
        velocity,
        magnetic,
        np.sin(math.pi * z),
        (1.0, 1.0, 2.0),
        math.pi,
    )
    assert transfer["direct_real_space"] > 0.0
    assert transfer["shell_sum"] == pytest.approx(transfer["direct_real_space"], abs=1.0e-14)
    assert transfer["closure_error"] == pytest.approx(0.0, abs=1.0e-14)
    assert transfer["transfer"][0] == pytest.approx(transfer["direct_real_space"])
    assert np.max(np.abs(transfer["transfer"][1:])) < 1.0e-14
    assert transfer["normalization_available"] is True

    context = products.BundleContext.__new__(products.BundleContext)
    x, y = products.curve_from_product(
        context,
        {"result": "pass", "pressure_transfer": transfer},
        "pressure_transfer.transfer_normalized_by_total",
    )
    assert x.tolist() == transfer["k_perp"]
    assert y.tolist() == transfer["transfer_normalized_by_total"]


def test_pressure_transfer_rejects_nonpositive_density_and_unavailable_normalization():
    """Pressure transfer fails closed on invalid density or a zero normalization."""

    fields = synthetic_fields((8, 8, 8))
    velocity = [fields["velx"], fields["vely"], fields["velz"]]
    magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
    density = np.array(fields["dens"], copy=True)
    density[0, 0, 0] = 0.0
    with pytest.raises(products.UnsupportedProduct, match="positive density"):
        products.pressure_transfer(
            density,
            velocity,
            magnetic,
            fields["p_perp"] - fields["eint"],
            (1.0, 1.0, 2.0),
            math.pi,
        )

    context = products.BundleContext.__new__(products.BundleContext)
    with pytest.raises(products.UnsupportedProduct, match="normalization is unavailable"):
        products.curve_from_product(
            context,
            {
                "result": "pass",
                "pressure_transfer": {
                    "normalization_available": False,
                    "k_perp": [0.0, 1.0],
                    "transfer_normalized_by_total": None,
                },
            },
            "pressure_transfer.transfer_normalized_by_total",
        )


def test_local_field_eddy_anisotropy_is_deterministic_and_count_weighted():
    """Fixed sampling yields replay-identical finite Figure 5(b) products."""

    assert (
        products.EDDY_SAMPLES,
        products.EDDY_BINS,
        products.EDDY_SEED,
        products.EDDY_ANGLE_DEGREES,
    ) == (2_000_000, 24, 731, 15.0)
    velocity, magnetic = synthetic_eddy_fields()
    first = products.local_field_eddy_anisotropy(
        velocity, magnetic, (1.0, 1.0, 2.0), 200_000, 10, 731
    )
    second = products.local_field_eddy_anisotropy(
        velocity, magnetic, (1.0, 1.0, 2.0), 200_000, 10, 731
    )
    assert first == second
    assert first["available"] is True
    assert first["velocity_perp"]["available"] is True
    assert first["magnetic_perp"]["available"] is True

    ensemble = products.mean_eddy_anisotropy([first, second])
    assert ensemble["samples_retained"] == 2 * first["samples_retained"]
    for name in ("velocity_perp", "magnetic_perp"):
        assert ensemble[name]["ell_perp_over_lperp"] == first[name]["ell_perp_over_lperp"]
        assert ensemble[name]["ell_parallel_over_lperp"] == pytest.approx(
            first[name]["ell_parallel_over_lperp"]
        )
        assert ensemble[name]["parallel_sample_counts"] == pytest.approx(
            2.0 * np.asarray(first[name]["parallel_sample_counts"])
        )
    context = products.BundleContext.__new__(products.BundleContext)
    x, y = products.curve_from_product(
        context,
        {"result": "pass", "eddy_anisotropy": ensemble},
        "eddy_anisotropy.velocity_perp",
    )
    assert x.tolist() == ensemble["velocity_perp"]["ell_perp_over_lperp"]
    assert y.tolist() == ensemble["velocity_perp"]["ell_parallel_over_lperp"]


def test_eddy_anisotropy_rejects_nonfinite_input_and_unavailable_curve():
    """Hostile fields and unavailable equal-power inversions fail closed."""

    velocity, magnetic = synthetic_eddy_fields((8, 8, 8))
    velocity[0] = np.array(velocity[0], copy=True)
    velocity[0][0, 0, 0] = np.nan
    with pytest.raises(products.UnsupportedProduct, match="nonfinite"):
        products.local_field_eddy_anisotropy(
            velocity, magnetic, (1.0, 1.0, 2.0), 10_000, 8, 731
        )

    context = products.BundleContext.__new__(products.BundleContext)
    with pytest.raises(products.UnsupportedProduct, match="do not overlap"):
        products.curve_from_product(
            context,
            {
                "result": "pass",
                "eddy_anisotropy": {
                    "velocity_perp": {
                        "available": False,
                        "reason": "structure functions do not overlap",
                    }
                },
            },
            "eddy_anisotropy.velocity_perp",
        )


def test_new_product_reference_semantics_are_stage_i_provenance_bound(
    tmp_path, monkeypatch
):
    """A product-name substitution cannot reuse an admitted reference identity."""

    references = tmp_path / "references"
    curve = references / "transfer.csv"
    write_curve(curve, [(0.25, 0.025, 0.01), (1.75, 0.175, 0.01)])
    manifest = references / "curves.json"
    write_json(manifest, {
        "schema_version": 1,
        "curves": [
            {
                "id": "fixture_transfer",
                "case": "fixture_case",
                "product": "pressure_transfer.transfer_normalized_by_total",
                "data_file": curve.name,
                "data_sha256": sha256(curve),
                "interpolation": "linear",
            }
        ],
        "surfaces": [],
    })
    stage_manifest = {
        "panel_status": {
            "analysis_case_aliases": {},
            "reference_manifests": {
                "fixture": {"path": manifest.name, "sha256": sha256(manifest)}
            },
            "reference_product_bindings": {
                "fixture_transfer": {
                    "kind": "curve",
                    "case": "fixture_case",
                    "product": "pressure_transfer.transfer_normalized_by_total",
                    "data_file": curve.name,
                    "data_sha256": sha256(curve),
                    "reference_manifest": "fixture",
                }
            },
            "panels": [{"id": "fig7bottom", "reference_products": ["fixture_transfer"]}],
        }
    }
    context = products.BundleContext.__new__(products.BundleContext)
    context.stage_manifest = stage_manifest
    context.reference_root = references.resolve()
    context.case_name = "fixture_case"
    ensemble = {
        "result": "pass",
        "pressure_transfer": {
            "normalization_available": True,
            "k_perp": [0.0, 1.0, 2.0],
            "transfer_normalized_by_total": [0.0, 0.1, 0.2],
        },
    }
    comparisons, bindings, deferred = products.reference_products_for_case(context, ensemble)
    comparison = comparisons["fixture_transfer"]
    assert comparison["available"] is True
    assert comparison["stage_i_binding_validated"] is True
    assert comparison["reference_data_file"] == curve.name
    assert comparison["reference_manifest_sha256"] == sha256(manifest)
    assert comparison["data_file"] == str(curve.resolve())
    assert comparison["data_sha256"] == sha256(curve)
    assert comparison["reference_y"] == comparison["reference_values"]
    assert comparison["simulated_y"] == comparison["simulated_values"]
    assert comparison["rms_normalized_by_reported_uncertainty"] == pytest.approx(
        comparison["normalized_residual_rms"]
    )
    assert len(bindings) == 2
    assert deferred == []
    binding = stage_manifest["panel_status"]["reference_product_bindings"][
        "fixture_transfer"
    ]
    monkeypatch.setattr(
        acceptance,
        "verified_reference_product",
        lambda *_: {
            "manifest_sha256": sha256(manifest),
            "data_path": str(curve.resolve()),
            "interpolation": "linear",
            "x": [0.25, 1.75],
            "reference": [0.025, 0.175],
            "uncertainty": [0.01, 0.01],
        },
    )
    rms, maximum, simulated = acceptance.normalized_reference_metrics(
        {}, "fixture_transfer", comparison, binding, "fixture_case"
    )
    assert rms == pytest.approx(0.0)
    assert maximum == pytest.approx(0.0)
    assert simulated == pytest.approx([0.025, 0.175])

    stage_manifest["panel_status"]["reference_product_bindings"]["fixture_transfer"][
        "product"
    ] = "eddy_anisotropy.velocity_perp"
    with pytest.raises(products.ScientificProductsError, match="differs from Stage I binding"):
        products.reference_products_for_case(context, ensemble)


def test_snapshot_analysis_rejects_mixed_authenticated_grids(monkeypatch):
    """Snapshots with differing physical grids cannot be averaged into one product."""

    def group(time: float) -> products.SnapshotGroup:
        return products.SnapshotGroup(
            time=time,
            representative=Path(f"/retained/{time:g}/rank_00000000/snapshot.bin"),
            rank_files=(),
            segment=f"s{int(time)}",
        )

    context = products.BundleContext.__new__(products.BundleContext)
    context.snapshots = [group(8.0), group(9.0)]
    context.case_name = "fixture_case"
    context.stage_manifest = {
        "panel_status": {
            "analysis_case_aliases": {},
            "reference_product_bindings": {},
        }
    }
    fields = synthetic_fields((8, 8, 8))
    monkeypatch.setattr(
        products,
        "read_snapshot_group",
        lambda value, _: (
            {name: np.array(field, copy=True) for name, field in fields.items()},
            (1.0, 1.0, 2.0 if value.time == 8.0 else 3.0),
            [],
        ),
    )
    ensemble, bindings, declared = products.analyze_snapshots(
        context, 8.0, 10.0, 16, [1, 2], 1_000_000, "auto"
    )
    assert ensemble["result"] == "inconclusive"
    assert "do not share one uniform analysis grid" in ensemble["reason"]
    assert bindings == []
    assert len(declared) == 2


def test_snapshot_analysis_derives_products_for_admitted_case_roles(monkeypatch):
    """Authenticated admitted roles trigger both new deterministic reductions."""

    group = products.SnapshotGroup(
        time=9.0,
        representative=Path("/retained/9/rank_00000000/snapshot.bin"),
        rank_files=(),
        segment="s09",
    )
    context = products.BundleContext.__new__(products.BundleContext)
    context.snapshots = [group]
    context.case_name = "fixture_case"
    context.stage_manifest = {
        "panel_status": {
            "analysis_case_aliases": {},
            "reference_product_bindings": {
                "fixture_transfer": {
                    "kind": "curve",
                    "case": "fixture_case",
                    "product": "pressure_transfer.transfer_normalized_by_total",
                },
                "fixture_eddy": {
                    "kind": "curve",
                    "case": "fixture_case",
                    "product": "eddy_anisotropy.velocity_perp",
                },
            },
        }
    }
    fields = synthetic_fields()
    velocity, magnetic = synthetic_eddy_fields()
    for name, value in zip(("velx", "vely", "velz"), velocity):
        fields[name] = value
    for name, value in zip(("bcc1", "bcc2", "bcc3"), magnetic):
        fields[name] = value
    monkeypatch.setattr(
        products,
        "read_snapshot_group",
        lambda *_: (
            {name: np.array(field, copy=True) for name, field in fields.items()},
            (1.0, 1.0, 2.0),
            [],
        ),
    )
    monkeypatch.setattr(products, "EDDY_SAMPLES", 200_000)
    monkeypatch.setattr(products, "EDDY_BINS", 10)
    ensemble, bindings, declared = products.analyze_snapshots(
        context, 8.0, 10.0, 16, [1, 2, 4], 1_000_000, "auto"
    )
    assert ensemble["result"] == "pass"
    assert ensemble["pressure_transfer"]["normalization_available"] is True
    assert ensemble["eddy_anisotropy"]["velocity_perp"]["available"] is True
    assert ensemble["admitted_curve_products"] == [
        "eddy_anisotropy.velocity_perp",
        "pressure_transfer.transfer_normalized_by_total",
    ]
    assert bindings == []
    assert len(declared) == 1


def test_snapshot_analysis_reports_reducer_rejection_as_inconclusive(monkeypatch):
    """A hostile ensemble reduction cannot escape as an unstructured failure."""

    group = products.SnapshotGroup(
        time=9.0,
        representative=Path("/retained/9/rank_00000000/snapshot.bin"),
        rank_files=(),
        segment="s09",
    )
    context = products.BundleContext.__new__(products.BundleContext)
    context.snapshots = [group]
    context.case_name = "fixture_case"
    context.stage_manifest = {
        "panel_status": {
            "analysis_case_aliases": {},
            "reference_product_bindings": {
                "fixture_transfer": {
                    "kind": "curve",
                    "case": "fixture_case",
                    "product": "pressure_transfer.transfer_normalized_by_total",
                }
            },
        }
    }
    fields = synthetic_fields((8, 8, 8))
    monkeypatch.setattr(
        products,
        "read_snapshot_group",
        lambda *_: (
            {name: np.array(field, copy=True) for name, field in fields.items()},
            (1.0, 1.0, 2.0),
            [],
        ),
    )
    def reject_ensemble(_):
        raise products.UnsupportedProduct("hostile pressure-transfer ensemble")

    monkeypatch.setattr(products, "mean_pressure_transfer", reject_ensemble)
    ensemble, _, _ = products.analyze_snapshots(
        context, 8.0, 10.0, 16, [1, 2], 1_000_000, "auto"
    )
    assert ensemble["result"] == "inconclusive"
    assert ensemble["reason"] == "hostile pressure-transfer ensemble"


def test_snapshot_reductions_emit_physical_k_convergence():
    """Supported reductions produce physical-k curves over reviewed k_perp/pi."""

    fields = synthetic_fields()
    lengths = (1.0, 1.0, 2.0)
    dk = math.pi
    spectra = {
        "velocity": products.shell_spectrum(
            [fields["velx"], fields["vely"], fields["velz"]], lengths, dk
        ),
        "magnetic_fluctuation": products.shell_spectrum(
            [fields["bcc1"], fields["bcc2"], fields["bcc3"]], lengths, dk
        ),
    }
    alignment = products.alignment_histograms(fields, lengths, [1, 2, 4], 16)
    convergence, deferred = products.convergence_products({
        "result": "pass",
        "spectra": spectra,
        "alignment": alignment,
    })
    assert convergence["peak_alignment"]["x"] == pytest.approx(
        [2.0 * math.pi, 4.0 * math.pi]
    )
    assert convergence["velocity_spectrum_shape"]["x"][0] == pytest.approx(
        4.0 * math.pi
    )
    assert convergence["velocity_spectrum_shape"]["x"][-1] <= 24.0 * math.pi
    assert convergence["velocity_spectrum_shape"]["x_coordinate"] == "physical_k_perp"
    assert "magnetic_fluctuation_spectrum_shape" in convergence
    assert not deferred


def test_bundle_snapshot_target_allows_leaf_symlink_but_rejects_escape(
    tmp_path: Path,
):
    """Bundle inventory may link to retained snapshots, never through an escaped parent."""

    bundle = tmp_path / "bundle"
    retained = tmp_path / "retained/snapshot.bin"
    retained.parent.mkdir()
    retained.write_bytes(b"snapshot")
    bundle.mkdir()
    (bundle / "snapshot.bin").symlink_to(retained)
    assert products.bundle_snapshot_target(
        bundle, "snapshot.bin", "bundle snapshot"
    ) == retained.resolve()
    with pytest.raises(products.ScientificProductsError, match="contained bundle-relative"):
        products.bundle_snapshot_target(bundle, "../retained/snapshot.bin", "bundle snapshot")
    (bundle / "escaped").symlink_to(retained.parent, target_is_directory=True)
    with pytest.raises(products.ScientificProductsError, match="escapes its authenticated root"):
        products.bundle_snapshot_target(
            bundle, "escaped/snapshot.bin", "bundle snapshot"
        )


def test_distinct_snapshot_groups_match_bundle_time_deduplication():
    """Whole-case inventory retains the ordered first group at each physical time."""

    def group(time: float, path: str) -> products.SnapshotGroup:
        return products.SnapshotGroup(
            time=time,
            representative=Path(path),
            rank_files=(
                {"path": path, "size_bytes": 1, "sha256": "a" * 64},
            ),
            segment="fixture",
        )

    first = group(8.0, "/retained/a/rank_00000000/snapshot.bin")
    duplicate = group(8.0 + 5.0e-13, "/retained/b/rank_00000000/snapshot.bin")
    later = group(8.25, "/retained/c/rank_00000000/snapshot.bin")
    assert products.distinct_snapshot_groups([later, duplicate, first]) == [first, later]


def test_skip_preserves_declared_snapshot_inventory():
    """Skipping expensive products still reports every selected declaration."""

    group = products.SnapshotGroup(
        time=9.0,
        representative=Path("/retained/rank_00000000/snapshot.bin"),
        rank_files=(
            {
                "path": "/retained/rank_00000000/snapshot.bin",
                "size_bytes": 1,
                "sha256": "a" * 64,
            },
        ),
        segment="s01_fixture",
    )
    context = products.BundleContext.__new__(products.BundleContext)
    context.snapshots = [group]
    ensemble, bindings, declared = products.analyze_snapshots(
        context, 8.0, 10.0, 16, [8], 100, "skip"
    )
    assert ensemble["selected_snapshot_count"] == 1
    assert bindings == []
    assert declared == ensemble["declared_snapshot_inventory"]


def fake_raw(
    logical: tuple[int, int, int], value: float
) -> dict[str, object]:
    """Return one tiny rank-local binary-parser payload."""

    block_shape = (2, 2, 2)
    return {
        "header": ["fixture"],
        "time": 8.0,
        "cycle": 1,
        "var_names": list(products.REQUIRED_FIELDS),
        "Nx1": 4,
        "Nx2": 2,
        "Nx3": 2,
        "nvars": len(products.REQUIRED_FIELDS),
        "x1min": 0.0,
        "x1max": 1.0,
        "x2min": 0.0,
        "x2max": 1.0,
        "x3min": 0.0,
        "x3max": 1.0,
        "nx1_mb": 2,
        "nx2_mb": 2,
        "nx3_mb": 2,
        "nx1_out_mb": 2,
        "nx2_out_mb": 2,
        "nx3_out_mb": 2,
        "mb_logical": np.asarray([[*logical, 0]], dtype=np.int64),
        "mb_index": np.asarray([[0, 1, 0, 1, 0, 1]], dtype=np.int64),
        "mb_geometry": np.asarray([[0.0, 0.5, 0.0, 1.0, 0.0, 1.0]]),
        "mb_data": {
            name: [np.full(block_shape, value if name != "dens" else 1.0)]
            for name in products.REQUIRED_FIELDS
        },
        "n_mbs": 1,
    }


def test_exact_rank_local_reconstruction_and_mutation_rejection(tmp_path, monkeypatch):
    """The raw parser consumes only the exact authenticated contiguous rank set."""

    rank_files = []
    raws = {}
    for rank, logical in enumerate(((0, 0, 0), (1, 0, 0))):
        path = tmp_path / f"rank_{rank:08d}/snapshot.bin"
        path.parent.mkdir()
        path.write_bytes(f"rank-{rank}".encode())
        rank_files.append({
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
        raws[str(path.resolve())] = fake_raw(logical, float(rank + 1))
    monkeypatch.setattr(
        products.bin_convert,
        "read_binary",
        lambda path: deepcopy(raws[str(Path(path).resolve())]),
    )
    group = products.SnapshotGroup(
        time=8.0,
        representative=Path(rank_files[0]["path"]),
        rank_files=tuple(rank_files),
        segment="s00_fixture",
    )
    fields, lengths, bindings = products.read_snapshot_group(group, 100)
    assert fields["velx"].shape == (2, 2, 4)
    assert np.all(fields["velx"][:, :, :2] == 1.0)
    assert np.all(fields["velx"][:, :, 2:] == 2.0)
    assert lengths == (1.0, 1.0, 1.0)
    assert len(bindings) == 2

    Path(rank_files[1]["path"]).write_bytes(b"forged")
    with pytest.raises(products.ScientificProductsError, match="SHA-256 differs"):
        products.read_snapshot_group(group, 100)


def test_candidate_write_is_no_clobber(tmp_path):
    """Evidence candidates cannot overwrite an existing output."""

    value = products.seal_evidence({"schema_version": 1})
    path = tmp_path / "candidate.json"
    products.write_candidate(path, value)
    with pytest.raises(products.ScientificProductsError, match="already exists"):
        products.write_candidate(path, value)


def test_candidate_write_rejects_canonical_root():
    """The standalone generator never writes beneath canonical storage."""

    with pytest.raises(products.ScientificProductsError, match="canonical root"):
        products.write_candidate(
            products.CANONICAL_ROOT / "forbidden-scientific-products.json",
            {"result": "forbidden"},
        )
