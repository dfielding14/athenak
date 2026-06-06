"""Focused regressions for retained-snapshot mechanism diagnostics."""

import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
ANALYZER_PATH = REPO_ROOT / "scripts" / "analyze_cgl_lf_paper.py"


def load_module(name, path):
    """Import one repository script without depending on the test-suite cwd."""

    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def analyzer():
    return load_module("analyze_cgl_lf_paper_mechanism_test", ANALYZER_PATH)


def mechanism_fields(amplitude=1.0):
    """Return a periodic field with nonzero mechanism and transfer diagnostics."""

    lengths = (1.0, 1.0, 2.0)
    shape = (16, 8, 8)
    z = (np.arange(shape[0]) + 0.5) * lengths[2] / shape[0]
    y = (np.arange(shape[1]) + 0.5) * lengths[1] / shape[1]
    x = (np.arange(shape[2]) + 0.5) * lengths[0] / shape[2]
    zz, _yy, xx = np.meshgrid(z, y, x, indexing="ij")
    fields = {
        "dens": np.ones(shape),
        "velx": amplitude * np.sin(2.0 * math.pi * xx),
        "vely": np.zeros(shape),
        "velz": amplitude * np.sin(math.pi * zz),
        "eint": np.ones(shape),
        "p_perp": 1.0 + 0.1 * np.cos(math.pi * zz),
        "bcc1": np.zeros(shape),
        "bcc2": np.zeros(shape),
        "bcc3": np.ones(shape),
    }
    return fields, lengths


def common_joint_ranges():
    """Return shared ranges spanning the focused multi-snapshot fixture."""

    field_ranges = {
        "delta_p": (-0.2, 0.2),
        "b_grad_delta_p": (-1.0, 1.0),
        "bb_grad_velocity": (-20.0, 20.0),
        "dln_b_dt": (-30.0, 30.0),
        "signed_pressure_stress_power_density": (-3.0, 3.0),
    }
    ranges = {
        "parallel": ((-1.0, 1.0), (-1.0, 1.0)),
        "perpendicular": ((-1.0, 1.0), (-1.0, 1.0)),
    }
    for name, (x_field, y_field) in {
        "b_grad_delta_p_vs_delta_p": ("delta_p", "b_grad_delta_p"),
        "bb_grad_velocity_vs_delta_p": ("delta_p", "bb_grad_velocity"),
        "bb_grad_velocity_vs_b_grad_delta_p": (
            "b_grad_delta_p",
            "bb_grad_velocity",
        ),
        "dln_b_dt_vs_bb_grad_velocity": ("bb_grad_velocity", "dln_b_dt"),
        "signed_pressure_stress_power_density_vs_bb_grad_velocity": (
            "bb_grad_velocity",
            "signed_pressure_stress_power_density",
        ),
        "signed_pressure_stress_power_density_vs_delta_p": (
            "delta_p",
            "signed_pressure_stress_power_density",
        ),
    }.items():
        ranges[name] = (field_ranges[x_field], field_ranges[y_field])
    return ranges


def test_mechanism_fields_and_paths_are_explicit(analyzer):
    fields, lengths = mechanism_fields()
    velocity = [fields["velx"], fields["vely"], fields["velz"]]
    magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
    products = analyzer.velocity_gradient_products(velocity, magnetic, lengths)

    assert np.array_equal(
        products["dln_b_dt"],
        products["bb_grad_velocity"] - products["div_velocity"],
    )
    assert np.allclose(
        products["bb_grad_velocity"],
        analyzer.periodic_gradient(fields["velz"], lengths)[2],
    )
    scalar_fields = analyzer.pdf_fields(fields, lengths)
    assert np.array_equal(
        scalar_fields["signed_pressure_stress_power_density"],
        -scalar_fields["delta_p"] * scalar_fields["bb_grad_velocity"],
    )

    record = analyzer.analyze_fields(fields, lengths, 1.25, 16, [1])
    for name in analyzer.MECHANISM_FIELD_DEFINITIONS:
        assert name in record["pdf"]
        assert name in record["spectra"]
        assert record["spectra"][name]["field_definition"] == (
            analyzer.MECHANISM_FIELD_DEFINITIONS[name]
        )
        assert record["mechanism_diagnostics"]["fields"][name][
            "scale_resolved_path"
        ] == f"spectra.{name}"
    assert record["spectra"]["grad_parallel_delta_p"] == (
        record["spectra"]["b_grad_delta_p"]
    )
    assert max(record["spectra"]["dln_b_dt"]["power_per_dk"]) > 0.0
    assert max(record["spectra"]["b_grad_delta_p"]["power_per_dk"]) > 0.0
    joints = record["mechanism_joint_diagnostics"]
    assert set(joints["products"]) == set(analyzer.MECHANISM_JOINT_FIELDS)
    for name, (x_field, y_field) in analyzer.MECHANISM_JOINT_FIELDS.items():
        product = joints["products"][name]
        assert product["x_field"] == x_field
        assert product["y_field"] == y_field
        assert product["joint_pdf"]["sample_count"] == fields["dens"].size
        assert product["joint_pdf"]["binned_sample_count"] == fields["dens"].size
        assert sum(product["conditional_y_given_x"]["sample_count"]) == (
            fields["dens"].size
        )
        assert "causal" not in joints["scope"].lower()
    assert record["mechanism_diagnostics"]["joint_conditional_products"][
        "dln_b_dt_vs_bb_grad_velocity"
    ]["path"] == (
        "mechanism_joint_diagnostics.products.dln_b_dt_vs_bb_grad_velocity"
    )


def test_pressure_stress_transfer_is_explicitly_signed(analyzer):
    fields, lengths = mechanism_fields()
    record = analyzer.analyze_fields(fields, lengths, 0.0, 16, [1])
    transfer = record["pressure_transfer"]
    work = record["pressure_work_decomposition"]

    assert transfer["signed_transfer"] == transfer["transfer"]
    assert transfer["signed_transfer_normalized_by_total"] == (
        transfer["transfer_normalized_by_total"]
    )
    assert transfer["direct_real_space"] < 0.0
    assert transfer["shell_sum"] == pytest.approx(transfer["direct_real_space"])
    assert transfer["direct_real_space"] == pytest.approx(
        work["anisotropic_stress_power"]
    )
    assert "kinetic-energy gain" in transfer["sign_convention"]
    assert record["mechanism_diagnostics"]["signed_pressure_stress_transfer"][
        "curve_path"
    ] == "pressure_transfer.signed_transfer"

    passive_record = analyzer.analyze_fields(
        fields, lengths, 0.0, 16, [1], model_choices={"passive_delta": "true"}
    )
    passive = passive_record["pressure_transfer"]
    assert passive["applied_to_flow"] is False
    assert "not applied to flow evolution" in passive["interpretation"]
    assert passive_record["mechanism_joint_diagnostics"]["applied_to_flow"] is False


def test_ensemble_exposes_descriptive_snapshot_and_block_uncertainty(analyzer):
    records = {}
    for index, amplitude in enumerate((1.0, 2.0, 3.0, 4.0)):
        fields, lengths = mechanism_fields(amplitude)
        records[str(index)] = analyzer.analyze_fields(
            fields,
            lengths,
            float(index),
            16,
            [1],
            joint_ranges=common_joint_ranges(),
        )

    ensemble = analyzer.average_snapshot_records(records)
    assert ensemble["snapshot_times"] == [0.0, 1.0, 2.0, 3.0]
    spectrum_uncertainty = ensemble["spectra"]["dln_b_dt"]["uncertainty"]
    transfer_uncertainty = ensemble["pressure_transfer"]["uncertainty"][
        "signed_transfer"
    ]
    work_uncertainty = ensemble["pressure_work_decomposition"]["uncertainty"][
        "anisotropic_stress_power"
    ]
    for uncertainty in (
        spectrum_uncertainty,
        transfer_uncertainty,
        work_uncertainty,
    ):
        assert uncertainty["available"] is True
        assert uncertainty["snapshot_count"] == 4
        assert uncertainty["contiguous_blocks"]["block_count"] == 2
        assert uncertainty["contiguous_blocks"]["block_snapshot_counts"] == [2, 2]
        assert "not realization-to-realization" in uncertainty["scope"]

    direct_values = np.asarray([
        records[str(index)]["pressure_transfer"]["direct_real_space"]
        for index in range(4)
    ])
    direct_uncertainty = ensemble["pressure_transfer"]["uncertainty"][
        "direct_real_space"
    ]
    assert direct_uncertainty["equal_snapshot_standard_error"] == pytest.approx(
        np.std(direct_values, ddof=1) / math.sqrt(4)
    )
    assert ensemble["mechanism_diagnostics"]["fields"]["dln_b_dt"][
        "scale_resolved_uncertainty_path"
    ] == "spectra.dln_b_dt.uncertainty"
    assert ensemble["mechanism_diagnostics"]["signed_pressure_stress_transfer"][
        "direct_real_space_path"
    ] == "pressure_transfer.direct_real_space_mean"
    joint_ensemble = ensemble["mechanism_joint_diagnostics"]
    assert joint_ensemble["snapshot_count"] == 4
    relationship = joint_ensemble["products"][
        "signed_pressure_stress_power_density_vs_bb_grad_velocity"
    ]
    assert relationship["joint_pdf"]["sample_count_sum"] == 4 * (
        records["0"]["shape_z_y_x"][0]
        * records["0"]["shape_z_y_x"][1]
        * records["0"]["shape_z_y_x"][2]
    )
    assert relationship["joint_pdf"]["uncertainty"]["snapshot_count"] == 4
    assert sum(relationship["conditional_y_given_x"]["sample_count"]) == (
        relationship["joint_pdf"]["binned_sample_count_sum"]
    )
    json.dumps(ensemble, allow_nan=False)


def test_conditional_profile_reports_histogram_quantiles_and_sign_fractions(analyzer):
    x_values = np.asarray([0.25, 0.25, 0.25, 0.25, 1.25, 1.25, 1.25, 1.25])
    y_values = np.asarray([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    legacy_density = analyzer.joint_pdf(
        x_values, y_values, 2, ((0.0, 2.0), (-2.0, 2.0))
    )
    counted_density = analyzer.joint_pdf(
        x_values,
        y_values,
        2,
        ((0.0, 2.0), (-2.0, 2.0)),
        include_counts=True,
    )
    assert np.array_equal(legacy_density["density"], counted_density["density"])
    assert "bin_counts" not in legacy_density
    assert counted_density["bin_counts"] == [[3, 1], [0, 4]]

    product = {
        "x_edges": [0.0, 1.0, 2.0],
        "y_edges": [-2.0, 0.0, 2.0],
        "density": [[3.0 / 16.0, 1.0 / 16.0], [0.0, 4.0 / 16.0]],
        "bin_counts": [[3, 1], [0, 4]],
    }
    profile = analyzer.conditional_profile_from_joint_pdf(product)

    assert profile["sample_count"] == [4, 4]
    assert profile["response_median"] == [-1.0, 1.0]
    assert profile["response_positive_fraction"] == pytest.approx([0.25, 1.0])
    assert profile["response_negative_fraction"] == pytest.approx([0.75, 0.0])
    assert "y-bin centers" in profile["bin_reconstruction"]
    assert "not an independent-sample inference" in profile["scope"]
    incompatible = {**product, "x_edges": [0.0, 0.5, 2.0]}
    with pytest.raises(ValueError, match="identical shared bin edges"):
        analyzer.mean_joint_distribution([product, incompatible], [0.0, 1.0])


def test_joint_coordinate_discovery_preserves_partial_legacy_paths(analyzer):
    delta_p = np.asarray([1.0, 2.0])
    b_grad_delta_p = np.asarray([3.0, 4.0])

    assert analyzer.mechanism_joint_coordinates({}) == {}
    coordinates = analyzer.mechanism_joint_coordinates({
        "delta_p": delta_p,
        "b_grad_delta_p": b_grad_delta_p,
    })
    assert set(coordinates) == {"b_grad_delta_p_vs_delta_p"}
    assert coordinates["b_grad_delta_p_vs_delta_p"] == (
        delta_p,
        b_grad_delta_p,
    )


def test_uncertainty_and_transfer_aliases_preserve_single_snapshot_and_legacy_use(
    analyzer,
):
    uncertainty = analyzer.descriptive_snapshot_block_uncertainty([[-1.0, 2.0]], [3.0])
    assert uncertainty["available"] is False
    assert uncertainty["equal_snapshot_standard_error"] is None
    assert uncertainty["contiguous_blocks"]["available"] is False

    unavailable_normalization = analyzer.mechanism_diagnostic_index(
        uncertainty_available=True,
        ensemble=True,
        normalized_transfer_available=False,
    )["signed_pressure_stress_transfer"]
    assert unavailable_normalization["normalized_curve_available"] is False
    assert "normalized_curve_uncertainty_path" not in unavailable_normalization

    legacy = {
        "pressure_transfer": {
            "k_perp": [0.0, 1.0],
            "transfer": [-1.0, 2.0],
            "normalization_available": True,
            "transfer_normalized_by_total": [-0.5, 1.0],
        }
    }
    x_old, y_old = analyzer.analyzed_product_curve(
        legacy, "pressure_transfer.transfer"
    )
    x_signed, y_signed = analyzer.analyzed_product_curve(
        legacy, "pressure_transfer.signed_transfer"
    )
    _, normalized_signed = analyzer.analyzed_product_curve(
        legacy, "pressure_transfer.signed_transfer_normalized_by_total"
    )
    assert np.array_equal(x_old, x_signed)
    assert np.array_equal(y_old, y_signed)
    assert np.array_equal(normalized_signed, [-0.5, 1.0])


def test_uncertainty_and_ensemble_averaging_preserve_legacy_record_compatibility(
    analyzer,
):
    uncertainty = analyzer.descriptive_snapshot_block_uncertainty(
        [[-1.0, 2.0], [1.0, 4.0]],
        [3.0, 3.0],
    )
    assert uncertainty["duplicate_snapshot_times_present"] is True
    assert uncertainty["unique_snapshot_time_count"] == 1

    records = {}
    for index, amplitude in enumerate((1.0, 2.0)):
        fields, lengths = mechanism_fields(amplitude)
        record = analyzer.analyze_fields(
            fields,
            lengths,
            float(index),
            16,
            [1],
            joint_ranges=common_joint_ranges(),
        )
        del record["mechanism_diagnostics"]
        del record["mechanism_joint_diagnostics"]
        records[str(index)] = record

    ensemble = analyzer.average_snapshot_records(records)
    assert "mechanism_diagnostics" not in ensemble
    assert "mechanism_joint_diagnostics" not in ensemble
    assert ensemble["pressure_transfer"]["signed_transfer"] == (
        ensemble["pressure_transfer"]["transfer"]
    )
