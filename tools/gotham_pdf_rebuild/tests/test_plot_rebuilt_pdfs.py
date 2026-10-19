"""Focused tests for the comprehensive rebuilt-PDF plotter."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


TOOLS_DIR = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "plot_rebuilt_pdfs", TOOLS_DIR / "plot_rebuilt_pdfs.py"
)
assert SPEC is not None and SPEC.loader is not None
PLOTTER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PLOTTER
SPEC.loader.exec_module(PLOTTER)


def _write_product(
    root: Path,
    product: str,
    *,
    weight: str = "mass",
    time_code: float = 10.5,
    sequence: str = "00001",
) -> None:
    product_dir = root / product
    product_dir.mkdir(parents=True)
    header = (
        """# small dense reconstructed PDF
format = dense
distribution = global_rebuild
ndim = 3
weight = mass
total_bins = 64
variable_1 = coord_r
nbin1 = 2
bin1_min = 1.0
bin1_max = 4.0
scale1 = log
stride1 = 16
variable_2 = coord_abscostheta
nbin2 = 2
bin2_min = 0.0
bin2_max = 1.0
scale2 = linear
stride2 = 4
variable_3 = temperature_kelvin
nbin3 = 2
bin3_min = 10.0
bin3_max = 1000.0
scale3 = log
stride3 = 1
bin_edges_1 = 1.0 2.0 4.0
bin_edges_2 = 0.0 0.5 1.0
bin_edges_3 = 10.0 100.0 1000.0
"""
    ).replace("weight = mass", f"weight = {weight}")
    product_dir.joinpath("gotham.header.pdf").write_text(
        header,
        encoding="utf-8",
    )
    values = np.zeros((4, 4, 4), dtype=np.float64)
    values[1:3, 1:3, 1:3] = np.arange(1.0, 9.0).reshape(2, 2, 2)
    values[0, 1, 1] = 4.0
    with product_dir.joinpath(f"gotham.{sequence}.pdf").open("wb") as stream:
        np.asarray([time_code], dtype=np.float64).tofile(stream)
        values.ravel(order="C").tofile(stream)
    root.joinpath("rebuild_manifest.json").write_text(
        json.dumps(
            {
                "products": [{"id": product}],
                "source_input_dir": "/archive/res_4pc/phase2/bin",
                "source_sequence": sequence,
                "output_number": sequence,
            }
        ),
        encoding="utf-8",
    )


def _write_four_dimensional_product(
    root: Path,
    product: str,
    *,
    physical_variables: tuple[str, str],
    weight_variable: str,
    sequence: str = "00001",
) -> None:
    product_dir = root / product
    product_dir.mkdir(parents=True)
    variables = (
        "coord_r",
        "coord_abscostheta",
        *physical_variables,
    )
    def physical_axis(variable: str) -> tuple[str, tuple[float, float, float]]:
        if "velocity" in variable:
            return "symlog", (-1000.0, 0.0, 1000.0)
        return "log", (10.0, 100.0, 1000.0)

    third_scale, third_edges = physical_axis(physical_variables[0])
    fourth_scale, fourth_edges = physical_axis(physical_variables[1])
    scales = ("log", "linear", third_scale, fourth_scale)
    edges = (
        (1.0, 4.0, 16.0),
        (0.0, 0.5, 1.0),
        third_edges,
        fourth_edges,
    )
    strides = (64, 16, 4, 1)
    lines = [
        "# small dense reconstructed PDF",
        "format = dense",
        "distribution = global_rebuild",
        "ndim = 4",
        (
            f"weight = {weight_variable}"
            if weight_variable in {"mass", "volume"}
            else "weight = variable"
        ),
        "total_bins = 256",
    ]
    if weight_variable not in {"mass", "volume"}:
        lines.insert(6, f"weight_variable = {weight_variable}")
    for index, (variable, scale, axis_edges, stride) in enumerate(
        zip(variables, scales, edges, strides), start=1
    ):
        lines.extend(
            (
                f"variable_{index} = {variable}",
                f"nbin{index} = 2",
                f"bin{index}_min = {axis_edges[0]}",
                f"bin{index}_max = {axis_edges[-1]}",
                f"scale{index} = {scale}",
                f"stride{index} = {stride}",
                f"bin_edges_{index} = {' '.join(str(value) for value in axis_edges)}",
            )
        )
        if scale == "symlog":
            lines.append(f"linthresh{index} = 5.0")
    product_dir.joinpath("gotham.header.pdf").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    values = np.zeros((4, 4, 4, 4), dtype=np.float64)
    values[1:3, 1:3, 1:3, 1:3] = np.arange(1.0, 17.0).reshape(2, 2, 2, 2)
    with product_dir.joinpath(f"gotham.{sequence}.pdf").open("wb") as stream:
        np.asarray([10.5], dtype=np.float64).tofile(stream)
        values.ravel(order="C").tofile(stream)
    root.joinpath("rebuild_manifest.json").write_text(
        json.dumps(
            {
                "products": [{"id": product}],
                "source_input_dir": "/archive/res_4pc/phase2/bin",
                "source_sequence": sequence,
                "output_number": sequence,
            }
        ),
        encoding="utf-8",
    )


def _product_metadata(output_dir: Path) -> dict:
    paths = list((output_dir / "metadata").glob("product.*.json"))
    assert len(paths) == 1
    return json.loads(paths[0].read_text())


def test_plotter_renders_scientific_heatmaps_and_conditional_profiles(
    tmp_path: Path,
) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_product(input_dir, "test_product")
    output_dir = tmp_path / "plots"

    result = PLOTTER.main(
        [
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--workers",
            "1",
            "--dpi",
            "40",
        ]
    )

    assert result == 0
    metadata = json.loads((output_dir / "metadata.json").read_text())
    assert metadata["product_count"] == 1
    assert metadata["plot_count"] == 3
    product = metadata["products"][0]
    assert product["plot_count"] == 3
    assert product["excluded_bins"]["excluded_absolute_fraction"] == 0.1
    assert product["overflow_by_axis"][0]["underflow_absolute_fraction"] == 0.1
    assert product["overflow_by_axis"][1]["underflow_absolute_fraction"] == 0.0
    assert product["snapshot_identity"]["time_title"] == "t =    50.0 [Myr]"
    assert (output_dir / "index.html").is_file()
    png_names = sorted(path.name for path in (output_dir / "png").glob("*.png"))
    assert len(png_names) == 3
    assert png_names == [
        "temperature_costheta_radius_shells_mass_weighted.res_4pc.t00001050.000Myr.png",
        "temperature_radius_costheta_mass_weighted.res_4pc.t00001050.000Myr.png",
        "temperature_radius_theta_cuts_mass_weighted.res_4pc.t00001050.000Myr.png",
    ]


def test_weight_presentation_treats_net_flux_as_signed() -> None:
    weight = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "edot_sph"}
    )

    assert weight.signed is True
    assert weight.style_key == "signed_net_energy_flux"
    assert weight.label == r"$\dot{E}_{\mathrm{net}}$"
    assert weight.factor > 0.0


def test_weight_presentation_keeps_cooling_luminosity_in_physical_units() -> None:
    weight = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "edot_cool"}
    )

    assert weight.signed is True
    assert weight.style_key == "signed_net_energy_flux"
    assert weight.label == r"$\dot{E}_{\mathrm{cool,net}}$"
    assert weight.unit == r"$\mathrm{erg\,s^{-1}}$"
    assert weight.factor == 1.0


def test_flux_labels_use_compact_physical_symbols() -> None:
    mdot = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "mdot_sph_out"}
    )
    edot = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "edot_sph_kin_in_abs"}
    )
    dust = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "dust_mdot_in_abs"}
    )

    assert mdot.label == r"$\dot{M}_{\mathrm{out}}$"
    assert edot.label == r"$|\dot{E}_{\mathrm{kin},\mathrm{in}}|$"
    assert dust.label == r"$|\dot{M}_{\mathrm{dust},\mathrm{in}}|$"


def test_radial_crossing_flux_labels_do_not_reintroduce_delta_r() -> None:
    mdot = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "mdot_sph_out"}
    )
    crossing = PLOTTER._crossing_presentation(mdot, "mdot_sph_out")

    assert crossing.label == r"$\dot{M}_{\mathrm{out}}$"
    assert crossing.unit == r"$M_\odot\,\mathrm{yr^{-1}}$"


def test_radial_ram_crossing_uses_momentum_flux_label() -> None:
    ram = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "radial_ram_out"}
    )
    crossing = PLOTTER._crossing_presentation(ram, "radial_ram_out")

    assert crossing.label == r"$\dot{p}_{\mathrm{out}}$"
    assert crossing.unit == r"$\mathrm{erg\,kpc^{-1}}$"


def test_thermokinematic_velocity_signs_follow_transport_direction() -> None:
    assert PLOTTER._thermokinematic_velocity_signs("edot_sph_out") == ("positive",)
    assert PLOTTER._thermokinematic_velocity_signs("edot_sph_in_abs") == (
        "negative",
    )
    assert PLOTTER._thermokinematic_velocity_signs("edot_sph") == (
        "positive",
        "negative",
    )


def test_fixed_radius_profiles_use_differential_axis_labels() -> None:
    mass = PLOTTER.weight_presentation({"weight": "mass"})
    temperature = PLOTTER.DisplayAxis(
        "temperature_kelvin",
        r"$T\;[\mathrm{K}]$",
        "log",
        np.asarray([1.0e3, 1.0e4, 1.0e5]),
        np.asarray([3.16227766e3, 3.16227766e4]),
        1.0,
    )
    costheta = PLOTTER.DisplayAxis(
        "coord_abscostheta",
        r"$|\cos\theta|$",
        "linear",
        np.asarray([0.0, 0.5, 1.0]),
        np.asarray([0.25, 0.75]),
        1.0,
    )

    assert PLOTTER._differential_profile_axis_label(
        mass, temperature
    ) == r"$\frac{\mathrm{d}M}{\mathrm{d}\log_{10} T}$"
    assert PLOTTER._differential_profile_axis_label(
        mass, costheta
    ) == r"$\frac{\mathrm{d}M}{\mathrm{d}|\cos\theta|}$"
    np.testing.assert_allclose(
        PLOTTER._differential_bin_widths(temperature), np.asarray([1.0, 1.0])
    )
    np.testing.assert_allclose(
        PLOTTER._differential_bin_widths(costheta), np.asarray([0.5, 0.5])
    )


def test_weight_presentation_distinguishes_inflow_and_enthalpy() -> None:
    inflow = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "mdot_sph_in_abs"}
    )
    enthalpy = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "edot_sph_th"}
    )

    assert inflow.style_key == "positive_inflow_mass_flux"
    assert inflow.signed is False
    assert r"\mathrm{enth}" in enthalpy.label
    assert inflow.filename_qualifier == "inflow_mass_flux"
    assert enthalpy.filename_qualifier == "net_enthalpy_flux"


def test_visual_norm_masks_zero_and_negligible_support() -> None:
    _, plotted, metadata = PLOTTER._weight_norm(
        np.asarray([0.0, 1.0, 1.0e6]), signed=True
    )

    assert plotted.mask.tolist() == [True, True, False]
    assert metadata["omitted_absolute_fraction"] < 1.0e-4


def test_movie_time_tokens_are_zero_padded_and_sort_chronologically() -> None:
    tokens = [PLOTTER._filename_time_token(value) for value in (9.5, 10.0, 10.5)]

    assert tokens == [
        "t00000950.000Myr",
        "t00001000.000Myr",
        "t00001050.000Myr",
    ]
    assert tokens == sorted(tokens)


def test_uniform_branch_has_distinct_descriptive_simulation_name(
    tmp_path: Path,
) -> None:
    input_dir = tmp_path / "res_8pc" / "phase2_uniform" / "00004"
    manifest = {
        "source_input_dir": "/archive/res_8pc/phase2_uniform/bin",
        "source_sequence": "00004",
        "output_number": "00004",
    }
    config = {
        "constants": {
            "feedback_start_myr": 1000.0,
            "time_title_precision_myr": 1,
            "time_title_numeric_width": 7,
        }
    }

    identity = PLOTTER._snapshot_identity(
        input_dir, manifest, time_code=10.0, config=config
    )

    assert identity["simulation"] == "res_8pc_uniform"
    assert identity["phase"] == "phase2_uniform"


def test_output24_profile_uses_short_semantic_movie_name() -> None:
    axes = [
        SimpleNamespace(raw_variable="coord_r"),
        SimpleNamespace(raw_variable="coord_abscostheta"),
        SimpleNamespace(raw_variable="hydro_w_e"),
    ]
    volume = PLOTTER.weight_presentation({"weight": "volume"})

    name = PLOTTER._semantic_plot_name(
        ["pressure", "radius", "costheta"], axes, {0, 1, 2}, volume
    )

    assert name == "pressure_radius_costheta"


def test_semantic_names_distinguish_signed_angle_and_transport_physics() -> None:
    folded_axes = [
        SimpleNamespace(raw_variable="coord_r"),
        SimpleNamespace(raw_variable="coord_abscostheta"),
    ]
    signed_axes = [
        SimpleNamespace(raw_variable="coord_r"),
        SimpleNamespace(raw_variable="coord_costheta"),
    ]
    outflow = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "edot_sph_out"}
    )
    kinetic = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "edot_sph_kin_out"}
    )

    folded = PLOTTER._semantic_plot_name(
        ["radius", "costheta"], folded_axes, {0, 1}, outflow
    )
    signed = PLOTTER._semantic_plot_name(
        ["radius", "signed_costheta"], signed_axes, {0, 1}, outflow
    )
    kinetic_name = PLOTTER._semantic_plot_name(
        ["radius", "costheta"], folded_axes, {0, 1}, kinetic
    )

    assert folded == "radius_costheta_outflow_energy_flux"
    assert signed == "radius_signed_costheta_outflow_energy_flux"
    assert kinetic_name == "radius_costheta_outflow_kinetic_energy_flux"


def test_semantic_names_record_marginalized_geometry() -> None:
    axes = [
        SimpleNamespace(raw_variable="coord_r"),
        SimpleNamespace(raw_variable="coord_abscostheta"),
        SimpleNamespace(raw_variable="temperature_kelvin"),
        SimpleNamespace(raw_variable="hydro_w_d"),
    ]
    mass = PLOTTER.weight_presentation({"weight": "mass"})

    name = PLOTTER._semantic_plot_name(
        ["temperature", "density"], axes, {2, 3}, mass
    )

    assert name == (
        "temperature_density_radius_integrated_costheta_integrated_mass_weighted"
    )


def test_weighted_mean_median_uses_conditional_histogram_weights() -> None:
    weights = np.asarray([[1.0, 3.0, 0.0], [0.0, 0.0, 0.0]])
    centers = np.asarray([10.0, 100.0, 1000.0])

    mean, median = PLOTTER._weighted_mean_median(weights, centers)

    np.testing.assert_allclose(mean[0], 77.5)
    np.testing.assert_allclose(median[0], 40.0)
    assert np.isnan(mean[1])
    assert np.isnan(median[1])


def test_conventional_pair_order_places_radius_before_angle() -> None:
    axes = [
        SimpleNamespace(raw_variable="coord_abscostheta"),
        SimpleNamespace(raw_variable="coord_r"),
    ]

    assert PLOTTER._ordered_pair(0, 1, axes) == (1, 0)


def test_conventional_pair_order_places_cylindrical_radius_before_height() -> None:
    axes = [
        SimpleNamespace(raw_variable="absolute_z"),
        SimpleNamespace(raw_variable="cylindrical_radius"),
    ]

    assert PLOTTER._ordered_pair(0, 1, axes) == (1, 0)


def test_vertical_geometry_semantic_name_matches_x_then_y_order() -> None:
    axes = [
        SimpleNamespace(raw_variable="absolute_z"),
        SimpleNamespace(raw_variable="cylindrical_radius"),
    ]
    ram = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "vertical_ram_out"}
    )
    x_index, y_index = PLOTTER._ordered_pair(0, 1, axes)

    name = PLOTTER._semantic_plot_name(
        [
            PLOTTER._axis_quantity_name(axes[x_index].raw_variable),
            PLOTTER._axis_quantity_name(axes[y_index].raw_variable),
        ],
        axes,
        {x_index, y_index},
        ram,
    )

    assert name == "cylindrical_radius_absolute_z_vertical_outflow_ram_flux"


def test_mass_and_volume_weighted_radial_angle_pairs_are_redundant() -> None:
    radius = SimpleNamespace(raw_variable="coord_r")
    angle = SimpleNamespace(raw_variable="coord_abscostheta")
    volume = PLOTTER.weight_presentation({"weight": "volume"})
    mass = PLOTTER.weight_presentation({"weight": "mass"})
    transport = PLOTTER.weight_presentation(
        {"weight": "variable", "weight_variable": "radial_mdot_out"}
    )

    assert PLOTTER._redundant_state_geometry_pair(radius, angle, volume)
    assert PLOTTER._redundant_state_geometry_pair(radius, angle, mass)
    assert not PLOTTER._redundant_state_geometry_pair(radius, angle, transport)


def test_angular_sector_weights_partition_folded_angle_bins() -> None:
    axis = SimpleNamespace(
        raw_variable="coord_abscostheta",
        edges=np.asarray([0.0, 0.5, 0.875, 1.0]),
    )
    equatorial = PLOTTER._angular_bin_selection_weights(axis, 0.0, 0.5)
    intermediate = PLOTTER._angular_bin_selection_weights(
        axis, 0.5, np.cos(np.deg2rad(30.0))
    )
    polar = PLOTTER._angular_bin_selection_weights(
        axis, np.cos(np.deg2rad(30.0)), 1.0
    )

    np.testing.assert_allclose(equatorial + intermediate + polar, np.ones(3))
    assert 0.0 < intermediate[1] < 1.0
    assert 0.0 < polar[1] < 1.0


def test_angular_sector_weights_fold_signed_costheta_hemispheres() -> None:
    axis = SimpleNamespace(
        raw_variable="coord_costheta",
        edges=np.asarray([-1.0, -0.5, 0.0, 0.5, 1.0]),
    )

    polar = PLOTTER._angular_bin_selection_weights(axis, 0.5, 1.0)

    np.testing.assert_allclose(polar, [1.0, 0.0, 0.0, 1.0])


def test_radial_bin_width_normalization() -> None:
    axis = SimpleNamespace(edges=np.asarray([1.0, 2.0, 4.0]))
    values = np.asarray([[2.0, 4.0], [8.0, 16.0]])

    normalized = PLOTTER._radial_bin_width_normalize(
        values, [axis, SimpleNamespace()], radial_axis_index=0
    )

    np.testing.assert_allclose(normalized, [[2.0, 4.0], [4.0, 8.0]])


def test_spherical_shell_volume_normalization() -> None:
    axis = SimpleNamespace(edges=np.asarray([1.0, 2.0, 4.0]))
    shell_volumes = 4.0 * np.pi / 3.0 * np.asarray([2.0**3 - 1.0, 4.0**3 - 2.0**3])
    values = np.asarray([[shell_volumes[0], 0.5 * shell_volumes[0]], [shell_volumes[1], 0.5 * shell_volumes[1]]])

    normalized = PLOTTER._spherical_shell_volume_normalize(
        values, [axis, SimpleNamespace()], radial_axis_index=0
    )

    np.testing.assert_allclose(normalized, [[1.0, 0.5], [1.0, 0.5]])


def test_radial_angular_wedge_volume_normalization() -> None:
    radial = SimpleNamespace(edges=np.asarray([1.0, 2.0]))
    angle = SimpleNamespace(
        raw_variable="coord_abscostheta",
        edges=np.asarray([0.0, 0.25, 1.0]),
    )
    shell_volume = 4.0 * np.pi / 3.0 * (2.0**3 - 1.0)
    values = np.asarray(
        [
            [0.25 * shell_volume, 0.125 * shell_volume],
            [0.75 * shell_volume, 0.375 * shell_volume],
        ]
    )

    normalized = PLOTTER._radial_angular_wedge_volume_normalize(
        values,
        angle,
        radial,
        angular_axis_index=0,
        radial_bin_index=0,
    )

    np.testing.assert_allclose(normalized, [[1.0, 0.5], [1.0, 0.5]])


@pytest.mark.parametrize("style_key", sorted(PLOTTER.NONNEGATIVE_HEATMAP_STYLE_KEYS))
def test_nonnegative_heatmap_colormaps_have_white_low_end(style_key: str) -> None:
    cmap = PLOTTER._configured_cmap({}, style_key)

    np.testing.assert_allclose(cmap(0.0)[:3], [1.0, 1.0, 1.0])


def test_nonnegative_heatmap_reverses_dark_low_configured_cmap() -> None:
    cmap = PLOTTER._configured_cmap(
        {"heatmap_weights": {"mass": {"cmap": "cmr.rainforest"}}},
        "mass",
    )

    assert cmap.name == "cmr.rainforest_r"
    np.testing.assert_allclose(cmap(0.0)[:3], [1.0, 1.0, 1.0])


def test_signed_heatmap_preserves_diverging_configured_cmap() -> None:
    cmap = PLOTTER._configured_cmap(
        {"heatmap_weights": {"signed_net_mass_flux": {"cmap": "RdBu_r"}}},
        "signed_net_mass_flux",
    )

    assert cmap.name == "RdBu_r"
    assert not np.allclose(cmap(0.0)[:3], [1.0, 1.0, 1.0])


def test_symlog_heatmap_colormap_uses_standard_white_centered_map() -> None:
    cmap = PLOTTER._heatmap_cmap(
        {},
        "signed_net_mass_flux",
        PLOTTER.colors.SymLogNorm(linthresh=1.0, vmin=-100.0, vmax=100.0),
    )

    np.testing.assert_allclose(cmap(0.5)[:3], [1.0, 1.0, 1.0])
    assert cmap.name == "cmr.fusion"


def test_volume_radial_panels_use_shell_volume_fraction(tmp_path: Path) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_product(input_dir, "test_product", weight="volume")
    output_dir = tmp_path / "plots"

    assert (
        PLOTTER.main(
            [
                str(input_dir),
                "--output-dir",
                str(output_dir),
                "--workers",
                "1",
                "--dpi",
                "40",
            ]
        )
        == 0
    )

    metadata = _product_metadata(output_dir)
    profile = next(
        plot
        for plot in metadata["plots"]
        if plot["kind"] == "conditional_mean_median_profiles"
    )
    assert profile["physical_axis"] == "temperature_kelvin"
    assert profile["statistics"] == ["mean", "median"]
    assert [series["sector"] for series in profile["radial_series"]] == [
        "all",
        "polar30",
        "equatorial30",
        "midlat",
    ]
    assert [
        series["target_radius_kpc"] for series in profile["angular_series"]
    ] == [2.0, 4.0, 8.0, 16.0, 32.0, 64.0]
    assert metadata["skipped_heatmaps"] == [
        {
            "reason": (
                "Mass- and volume-weighted radius--polar-angle marginals "
                "reproduce weighted geometry rather than the physical "
                "quantity carried by the product."
            ),
            "x_axis": "coord_r",
            "y_axis": "coord_abscostheta",
        }
    ]
    radial_heatmap = next(
        plot
        for plot in metadata["plots"]
        if plot["kind"] == "angular_sector_radial_marginals"
    )
    shell_heatmap = next(
        plot
        for plot in metadata["plots"]
        if plot["kind"] == "fixed_radial_shell_angular_marginals"
    )
    assert radial_heatmap["radial_normalization"] == (
        "divided_by_angular_sector_shell_volume"
    )
    assert shell_heatmap["radial_normalization"] == (
        "divided_by_radial_angular_wedge_volume"
    )
    assert [panel["sector"] for panel in radial_heatmap["panels"]] == [
        "all",
        "polar30",
        "equatorial30",
        "midlat",
    ]
    assert [panel["target_radius_kpc"] for panel in shell_heatmap["panels"]] == [
        2.0,
        4.0,
        8.0,
        16.0,
        32.0,
        64.0,
    ]
    assert radial_heatmap["weight_unit"] == ""
    assert shell_heatmap["weight_unit"] == ""
    assert radial_heatmap["cmap"].endswith("_r")
    assert shell_heatmap["cmap"].endswith("_r")
    assert metadata["weight"]["radial_volume_note"]


def test_conditioned_phase_product_makes_only_four_shell_figures(
    tmp_path: Path,
) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_four_dimensional_product(
        input_dir,
        "science_phase_density_temperature_volume",
        physical_variables=(
            "hydrogen_number_density_cm3",
            "temperature_kelvin",
        ),
        weight_variable="volume",
    )
    output_dir = tmp_path / "plots"

    assert (
        PLOTTER.main(
            [
                str(input_dir),
                "--output-dir",
                str(output_dir),
                "--workers",
                "1",
                "--dpi",
                "40",
            ]
        )
        == 0
    )

    metadata = _product_metadata(output_dir)
    assert metadata["plot_count"] == 4
    assert {
        plot["kind"] for plot in metadata["plots"]
    } == {"conditioned_physical_pair_radial_shells"}
    png_names = [path.name for path in (output_dir / "png").glob("*.png")]
    assert len(png_names) == 4
    assert all("output" not in name for name in png_names)
    assert not any("radius_costheta" in name for name in png_names)


def test_thermokinematic_product_makes_split_velocity_dedicated_figures(
    tmp_path: Path,
) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_four_dimensional_product(
        input_dir,
        "science_transport_thermokinematic_mdot_out",
        physical_variables=("temperature_kelvin", "velocity_r_km_s"),
        weight_variable="mdot_sph_out",
    )
    output_dir = tmp_path / "plots"

    assert (
        PLOTTER.main(
            [
                str(input_dir),
                "--output-dir",
                str(output_dir),
                "--workers",
                "1",
                "--dpi",
                "40",
            ]
        )
        == 0
    )

    metadata = _product_metadata(output_dir)
    assert metadata["plot_count"] == 7
    assert sum(
        plot["kind"] == "thermokinematic_velocity_temperature_radial_shells"
        for plot in metadata["plots"]
    ) == 4
    png_names = [path.name for path in (output_dir / "png").glob("*.png")]
    assert len(png_names) == 7
    assert all("output" not in name for name in png_names)
    assert sum("positive_vr" in name for name in png_names) == 4
    assert sum("negative_vr" in name for name in png_names) == 0
    assert all(
        "radial_velocity_temperature_radius_shells" in name
        for name in png_names
        if "positive_vr" in name or "negative_vr" in name
    )


def test_source_context_marks_unvalidated_reduction() -> None:
    context = PLOTTER._source_context({"shards_processed": 8})

    assert context["status"] == "partial_or_unvalidated"
    assert context["annotation"] == "PARTIAL / UNVALIDATED REDUCTION: 8 shards"


def test_overwrite_replaces_only_matching_product_outputs(tmp_path: Path) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_product(input_dir, "test_product")
    output_dir = tmp_path / "plots"
    stale = output_dir / "png" / "unrelated.png"
    stale.parent.mkdir(parents=True)
    stale.write_text("obsolete", encoding="utf-8")
    common = [
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--workers",
        "1",
        "--dpi",
        "40",
    ]
    assert PLOTTER.main(common + ["--overwrite"]) == 0
    product_path = Path(_product_metadata(output_dir)["plots"][0]["paths"][0])
    product_path.write_text("obsolete", encoding="utf-8")

    assert PLOTTER.main(common + ["--overwrite"]) == 0

    assert stale.is_file()
    assert product_path.stat().st_size > len("obsolete")
    assert len(list((output_dir / "png").glob("*.png"))) == 4


def test_skip_existing_invalidates_when_dpi_changes(tmp_path: Path) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_product(input_dir, "test_product")
    output_dir = tmp_path / "plots"
    common = [
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--workers",
        "1",
    ]

    assert PLOTTER.main(common + ["--dpi", "40"]) == 0
    stale = output_dir / "png" / "unrelated.png"
    stale.write_text("obsolete", encoding="utf-8")
    assert PLOTTER.main(common + ["--dpi", "50", "--skip-existing"]) == 0

    metadata = json.loads((output_dir / "metadata.json").read_text())
    assert metadata["products"][0]["cached"] is False
    assert metadata["products"][0]["fingerprint"]["dpi"] == 50
    assert stale.is_file()


def test_distinct_snapshots_append_to_one_flat_movie_directory(tmp_path: Path) -> None:
    first_input = tmp_path / "first"
    first_input.mkdir()
    _write_product(first_input, "test_product", time_code=10.5, sequence="00001")
    second_input = tmp_path / "second"
    second_input.mkdir()
    _write_product(second_input, "test_product", time_code=10.6, sequence="00002")
    output_dir = tmp_path / "plots"
    common = ["--output-dir", str(output_dir), "--workers", "1", "--dpi", "40"]

    assert PLOTTER.main([str(first_input), *common]) == 0
    assert PLOTTER.main([str(second_input), *common]) == 0

    png_names = sorted(path.name for path in (output_dir / "png").glob("*.png"))
    assert len(png_names) == 6
    assert any(".t00001050.000Myr." in name for name in png_names)
    assert any(".t00001060.000Myr." in name for name in png_names)
    assert len(list((output_dir / "metadata").glob("product.*.json"))) == 2


def test_plotter_refuses_nested_output_tree(tmp_path: Path) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_product(input_dir, "test_product")

    with pytest.raises(SystemExit, match="nested input/output trees"):
        PLOTTER.main(
            [
                str(input_dir),
                "--output-dir",
                str(input_dir / "plots"),
            ]
        )


def test_plotter_refuses_missing_completed_manifest(tmp_path: Path) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_product(input_dir, "test_product")
    (input_dir / "rebuild_manifest.json").unlink()

    with pytest.raises(ValueError, match="Completed reducer manifest is missing"):
        PLOTTER.main([str(input_dir)])


def test_plotter_reports_nonfinite_histogram_as_failure(tmp_path: Path) -> None:
    input_dir = tmp_path / "rebuild"
    input_dir.mkdir()
    _write_product(input_dir, "test_product")
    payload = input_dir / "test_product" / "gotham.00001.pdf"
    with payload.open("r+b") as stream:
        stream.seek(8 + 21 * 8)
        np.asarray([np.nan], dtype=np.float64).tofile(stream)
    output_dir = tmp_path / "plots"

    result = PLOTTER.main(
        [
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--workers",
            "1",
            "--dpi",
            "40",
        ]
    )

    assert result == 1
    metadata = json.loads((output_dir / "metadata.json").read_text())
    assert metadata["failed_product_count"] == 1
    assert "non-finite histogram weights" in metadata["failures"][0]["error"]
