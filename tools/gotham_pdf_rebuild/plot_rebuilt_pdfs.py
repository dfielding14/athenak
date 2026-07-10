#!/usr/bin/env python3
"""Plot informative marginals of every reconstructed GOTHAM PDF product.

For each product this writes every informative unique pair of axes, summing
over all other axes. Mass- and volume-weighted radius--polar-angle pairs are
omitted because they reproduce weighted geometry rather than the physical
quantity carried by the product. Products containing radius, polar angle, and
another physical variable replace the separate radial and angle heatmaps with
one four-panel radial heatmap split into all-angle, polar, midplane, and
intermediate sectors, one six-panel angular heatmap at fixed radial shells,
and one two-panel conditional mean/median line summary.

The presentation follows the GOTHAM analysis slice/PDF preview conventions:
configured colormaps, matched-height colorbars, mathtext labels, four-sided
inward ticks, physical axis conversions, fixed-width time titles, and
under/overflow records in metadata.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import fcntl
import fnmatch
import hashlib
import html
import json
import math
import os
import re
import sys
import traceback
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.lines import Line2D

try:
    import cmasher  # noqa: F401  # Registers cmr.* maps when installed.
except ImportError:
    cmasher = None


DEFAULT_ANALYSIS_ROOT = Path("/lustre/orion/ast207/proj-shared/gotham/analysis")
DEFAULT_STYLE_CONFIG = DEFAULT_ANALYSIS_ROOT / "config/plot_defaults.yaml"
PAYLOAD_RE = re.compile(r"^gotham\.(\d{5})\.pdf$")

M_SUN_CGS = 1.98847e33
YEAR_S = 3.15576e7
MYR_CGS = 3.15576e13
TIME_CGS = 3.15576e15
Z_SOLAR = 0.02
EMPTY_SUPPORT_COLOR = "#f7f7f7"
VISUAL_OMITTED_ABSOLUTE_WEIGHT_FRACTION = 1.0e-4
FILENAME_TIME_INTEGER_WIDTH = 8
FILENAME_TIME_DECIMALS = 3
SIMULATION_RE = re.compile(r"^res_(?:4|8)pc(?:_highmdot|_lowmdot)?$")

ANGULAR_SECTOR_ORDER = ("all", "polar30", "equatorial30", "midlat")
FIXED_RADIAL_SHELL_TARGETS_KPC = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
DISPLAY_AXIS_LIMITS = {
    "hydro_w_s_00": (0.1, 10.0),
}
THERMOKINEMATIC_VELOCITY_MIN_KM_S = 1.0
THERMOKINEMATIC_VELOCITY_MAX_KM_S = 1.0e4
THERMOKINEMATIC_MACH_GUIDES = (10.0**-0.5, 1.0, 10.0**0.5)
THERMOKINEMATIC_BERNOULLI_GUIDES_KM_S = (
    10.0,
    10.0**1.5,
    100.0,
    10.0**2.5,
    1000.0,
    10.0**3.5,
)
THERMOKINEMATIC_GAMMA = 5.0 / 3.0
THERMOKINEMATIC_SOUND_SPEED_KM_S_PER_SQRT_K = 0.14896055403315162
ANGULAR_SECTOR_TITLES = {
    "all": r"All $\theta$",
    "polar30": r"Within $30^\circ$ of the $z$ axis",
    "equatorial30": r"Within $30^\circ$ of the midplane",
    "midlat": r"Intermediate $30^\circ$",
}

SIGNED_WEIGHT_VARIABLES = {
    "mdot_sph",
    "mdot_sph_in",
    "edot_sph",
    "edot_sph_in",
    "edot_sph_kin",
    "edot_sph_th",
    "edot_cool",
}

NEGATIVE_INFLOW_WEIGHT_VARIABLES = {"mdot_sph_in", "edot_sph_in"}

CMAP_FALLBACKS = {
    "mass": "cmr.rainforest_r",
    "volume": "cmr.chroma_r",
    "volume_fraction": "cmr.chroma_r",
    "signed_net_mass_flux": "cmr.fusion",
    "signed_net_energy_flux": "cmr.fusion",
    "negative_inflow_mass_flux": "cmr.fusion",
    "negative_inflow_energy_flux": "cmr.fusion",
    "symlog_diverging": "cmr.fusion",
    "positive_inflow_mass_flux": "cmr.arctic_r",
    "positive_inflow_energy_flux": "cmr.amethyst_r",
    "positive_outflow_mass_flux": "cmr.torch_r",
    "positive_outflow_ram_momentum_flux": "cmr.voltage_r",
    "positive_outflow_energy_flux": "cmr.sunburst_r",
    "code": "cmr.freeze_r",
}

NONNEGATIVE_HEATMAP_STYLE_KEYS = {
    "mass",
    "volume",
    "volume_fraction",
    "positive_inflow_mass_flux",
    "positive_inflow_energy_flux",
    "positive_outflow_mass_flux",
    "positive_outflow_ram_momentum_flux",
    "positive_outflow_energy_flux",
    "code",
}

WEIGHT_LABELS = {
    "mdot_sph": r"$\dot{M}_{\mathrm{net}}$",
    "mdot_sph_out": r"$\dot{M}_{\mathrm{out}}$",
    "mdot_sph_in": r"$\dot{M}_{\mathrm{in}}$",
    "mdot_sph_in_abs": r"$|\dot{M}_{\mathrm{in}}|$",
    "edot_sph": r"$\dot{E}_{\mathrm{net}}$",
    "edot_sph_out": r"$\dot{E}_{\mathrm{out}}$",
    "edot_sph_in": r"$\dot{E}_{\mathrm{in}}$",
    "edot_sph_in_abs": r"$|\dot{E}_{\mathrm{in}}|$",
    "edot_sph_kin": r"$\dot{E}_{\mathrm{kin},\mathrm{net}}$",
    "edot_sph_kin_out": r"$\dot{E}_{\mathrm{kin},\mathrm{out}}$",
    "edot_sph_kin_in_abs": r"$|\dot{E}_{\mathrm{kin},\mathrm{in}}|$",
    "edot_sph_th": r"$\dot{E}_{\mathrm{enth},\mathrm{net}}$",
    "edot_sph_th_out": r"$\dot{E}_{\mathrm{enth},\mathrm{out}}$",
    "edot_sph_th_in_abs": r"$|\dot{E}_{\mathrm{enth},\mathrm{in}}|$",
    "edot_cool": r"$\dot{E}_{\mathrm{cool,net}}$",
    "radial_ram_out": r"$\sum \rho v_{r,+}^2\,dV$",
    "vertical_mdot_out": r"$\dot{M}_{\perp,\mathrm{out}}$",
    "vertical_mdot_in_abs": r"$|\dot{M}_{\perp,\mathrm{in}}|$",
    "vertical_edot_out": r"$\dot{E}_{\perp,\mathrm{out}}$",
    "vertical_ram_out": r"$\sum \rho v_{\perp,+}^2\,dV$",
    "gas_metal_mdot_out": r"$\dot{M}_{Z_{\mathrm{gas}},\mathrm{out}}$",
    "gas_metal_mdot_in_abs": r"$|\dot{M}_{Z_{\mathrm{gas}},\mathrm{in}}|$",
    "total_metal_mdot_out": r"$\dot{M}_{Z_{\mathrm{total}},\mathrm{out}}$",
    "total_metal_mdot_in_abs": r"$|\dot{M}_{Z_{\mathrm{total}},\mathrm{in}}|$",
    "dust_mdot_out": r"$\dot{M}_{\mathrm{dust},\mathrm{out}}$",
    "dust_mdot_in_abs": r"$|\dot{M}_{\mathrm{dust},\mathrm{in}}|$",
}

WEIGHT_FILENAME_QUALIFIERS = {
    "total_metal_mass": "total_metal_mass_weighted",
    "dust_mass": "dust_mass_weighted",
    "mdot_sph": "net_mass_flux",
    "mdot_sph_out": "outflow_mass_flux",
    "mdot_sph_in": "negative_inflow_mass_flux",
    "mdot_sph_in_abs": "inflow_mass_flux",
    "edot_sph": "net_energy_flux",
    "edot_sph_out": "outflow_energy_flux",
    "edot_sph_in": "negative_inflow_energy_flux",
    "edot_sph_in_abs": "inflow_energy_flux",
    "edot_sph_kin": "net_kinetic_energy_flux",
    "edot_sph_kin_out": "outflow_kinetic_energy_flux",
    "edot_sph_kin_in_abs": "inflow_kinetic_energy_flux",
    "edot_sph_th": "net_enthalpy_flux",
    "edot_sph_th_out": "outflow_enthalpy_flux",
    "edot_sph_th_in_abs": "inflow_enthalpy_flux",
    "edot_cool": "net_cooling_luminosity",
    "radial_ram_out": "outflow_radial_ram_flux",
    "vertical_mdot_out": "vertical_outflow_mass_flux",
    "vertical_mdot_in_abs": "vertical_inflow_mass_flux",
    "vertical_edot_out": "vertical_outflow_energy_flux",
    "vertical_ram_out": "vertical_outflow_ram_flux",
    "gas_metal_mdot_out": "outflow_gas_metal_mass_flux",
    "gas_metal_mdot_in_abs": "inflow_gas_metal_mass_flux",
    "total_metal_mdot_out": "outflow_total_metal_mass_flux",
    "total_metal_mdot_in_abs": "inflow_total_metal_mass_flux",
    "dust_mdot_out": "outflow_dust_mass_flux",
    "dust_mdot_in_abs": "inflow_dust_mass_flux",
}


@dataclass(frozen=True)
class DisplayAxis:
    """One PDF axis converted into display units."""

    raw_variable: str
    label: str
    scale: str
    edges: np.ndarray
    centers: np.ndarray
    linthresh: float


@dataclass(frozen=True)
class WeightPresentation:
    """Display semantics for one integrated histogram weight."""

    style_key: str
    label: str
    unit: str
    factor: float
    signed: bool
    filename_qualifier: Optional[str] = None


def _strip_matplotlib_prefix(name: str) -> str:
    prefix = "matplotlib:"
    return name[len(prefix) :] if name.startswith(prefix) else name


def _safe_component(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _json_ready(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_ready(item) for item in value]
    return value


def _write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(_json_ready(value), stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _acquire_output_lock(output_dir: Path) -> Any:
    lock_path = output_dir / ".plot_rebuilt_pdfs.lock"
    stream = lock_path.open("a+", encoding="utf-8")
    try:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        stream.close()
        raise SystemExit(f"Another plotter process holds the output lock: {lock_path}")
    stream.seek(0)
    stream.truncate()
    stream.write(f"pid={os.getpid()} host={os.uname().nodename}\n")
    stream.flush()
    return stream


@lru_cache(maxsize=4)
def _analysis_api(analysis_root_string: str) -> Tuple[Any, Any, Any]:
    analysis_root = Path(analysis_root_string)
    if not analysis_root.is_dir():
        raise FileNotFoundError(f"GOTHAM analysis root does not exist: {analysis_root}")
    if str(analysis_root) not in sys.path:
        sys.path.insert(0, str(analysis_root))
    from gotham_analysis.plotting.style import add_matched_colorbar, apply_axes_style
    from gotham_analysis.readers.pdf import read_pdf

    return read_pdf, add_matched_colorbar, apply_axes_style


def _time_title(time_code: float, config: Mapping[str, Any]) -> Tuple[str, float]:
    constants = config["constants"]
    time_since_feedback_myr = (
        float(time_code) * TIME_CGS / MYR_CGS
        - float(constants["feedback_start_myr"])
    )
    precision = int(constants["time_title_precision_myr"])
    rounded = round(time_since_feedback_myr, precision)
    if rounded == 0:
        rounded = 0.0
    title = (
        f"t = {rounded:{int(constants['time_title_numeric_width'])}.{precision}f} "
        "[Myr]"
    )
    return title, time_since_feedback_myr


def _filename_time_token(time_code: float) -> str:
    absolute_time_myr = float(time_code) * TIME_CGS / MYR_CGS
    if not math.isfinite(absolute_time_myr) or absolute_time_myr < 0.0:
        raise ValueError(
            "Movie-ready filename tokens require finite non-negative absolute time"
        )
    width = FILENAME_TIME_INTEGER_WIDTH + 1 + FILENAME_TIME_DECIMALS
    maximum = 10**FILENAME_TIME_INTEGER_WIDTH
    if absolute_time_myr >= maximum:
        raise ValueError(
            f"Absolute time {absolute_time_myr} Myr exceeds filename width"
        )
    return f"t{absolute_time_myr:0{width}.{FILENAME_TIME_DECIMALS}f}Myr"


def _first_matching_path_component(paths: Iterable[Path], pattern: Any) -> Optional[str]:
    for path in paths:
        for component in path.parts:
            if pattern(component):
                return component
    return None


def _snapshot_identity(
    input_dir: Path,
    manifest: Mapping[str, Any],
    *,
    time_code: float,
    config: Mapping[str, Any],
) -> Dict[str, Any]:
    source_input = Path(str(manifest.get("source_input_dir", "")))
    simulation = _first_matching_path_component(
        (input_dir, source_input), lambda value: SIMULATION_RE.match(value) is not None
    )
    phase = _first_matching_path_component(
        (input_dir, source_input), lambda value: value.startswith("phase")
    )
    if phase == "phase2_uniform":
        simulation = f"{simulation}_uniform" if simulation else "unknown_simulation_uniform"
    sequence = str(manifest.get("source_sequence") or input_dir.name)
    output_number = str(manifest.get("output_number") or sequence)
    title, time_since_feedback_myr = _time_title(time_code, config)
    return {
        "simulation": simulation or "unknown_simulation",
        "phase": phase or "unknown_phase",
        "source_sequence": sequence,
        "output_number": output_number,
        "time_code": float(time_code),
        "absolute_time_myr": float(time_code) * TIME_CGS / MYR_CGS,
        "time_since_feedback_myr": time_since_feedback_myr,
        "time_title": title,
        "time_token": _filename_time_token(time_code),
    }


def _identity_components(identity: Mapping[str, Any]) -> List[str]:
    return [
        _safe_component(str(identity["simulation"])),
        _safe_component(str(identity["time_token"])),
        _safe_component(str(identity["phase"])),
        "s" + _safe_component(str(identity["source_sequence"])),
        "o" + _safe_component(str(identity["output_number"])),
    ]


def _plot_stem(name: str, identity: Mapping[str, Any]) -> str:
    return ".".join(
        (
            _safe_component(name),
            _safe_component(str(identity["simulation"])),
            _safe_component(str(identity["time_token"])),
        )
    )


def _metadata_stem(
    kind: str, identity: Mapping[str, Any], *, product: Optional[str] = None
) -> str:
    components = [_safe_component(kind), *_identity_components(identity)[:2]]
    if product is not None:
        components.append(_safe_component(product))
    components.extend(_identity_components(identity)[2:])
    return ".".join(components)


def load_style_config(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a YAML mapping")
    return value


def _axis_factor_and_label(variable: str) -> Tuple[float, str]:
    labels = {
        "coord_r": (1.0, r"$r\;[\mathrm{kpc}]$"),
        "coord_abscostheta": (1.0, r"$|\cos\theta|$"),
        "coord_costheta": (1.0, r"$\cos\theta$"),
        "temperature_kelvin": (1.0, r"$T\;[\mathrm{K}]$"),
        "hydrogen_number_density_cm3": (1.0, r"$n_{\mathrm{H}}\;[\mathrm{cm^{-3}}]$"),
        "pressure_over_kb_k_cm3": (1.0, r"$P/k_B\;[\mathrm{K\,cm^{-3}}]$"),
        "cooling_rate_erg_s_cm3": (
            1.0,
            r"$\dot{e}_{\mathrm{cool,net}}\;[\mathrm{erg\,s^{-1}\,cm^{-3}}]$",
        ),
        "cooling_time_myr": (1.0, r"$t_{\mathrm{cool,net}}\;[\mathrm{Myr}]$"),
        "velocity_r_km_s": (1.0, r"$v_r\;[\mathrm{km\,s^{-1}}]$"),
        "velocity_theta_km_s": (1.0, r"$v_\theta\;[\mathrm{km\,s^{-1}}]$"),
        "velocity_phi_km_s": (1.0, r"$v_\phi\;[\mathrm{km\,s^{-1}}]$"),
        "absolute_velocity_r_km_s": (1.0, r"$|v_r|\;[\mathrm{km\,s^{-1}}]$"),
        "sound_speed_km_s": (1.0, r"$c_s\;[\mathrm{km\,s^{-1}}]$"),
        "mach": (1.0, r"$\mathcal{M}$"),
        "radial_mach": (1.0, r"$\mathcal{M}_r$"),
        "absolute_radial_mach": (1.0, r"$|\mathcal{M}_r|$"),
        "absolute_z": (1.0, r"$|z|\;[\mathrm{kpc}]$"),
        "cylindrical_radius": (1.0, r"$R_{\mathrm{cyl}}\;[\mathrm{kpc}]$"),
        "entropy_proxy": (1.0, r"$K_{\mathrm{proxy}}\;[\mathrm{code}]$"),
        "specific_angular_momentum_z": (1.0, r"$j_z\;[\mathrm{code}]$"),
        "specific_angular_momentum_z_kpc_km_s": (
            1.0,
            r"$j_z\;[\mathrm{kpc\,km\,s^{-1}}]$",
        ),
        "hydro_w_s_00": (1.0 / Z_SOLAR, r"$Z_{\mathrm{gas}}/Z_\odot$"),
        "hydro_w_s_01": (1.0, r"$D_{\mathrm{small}}$"),
        "hydro_w_s_02": (1.0, r"$D_{\mathrm{large}}$"),
        "dust_to_total_metal": (1.0, r"$D/Z_{\mathrm{total}}$"),
        "small_grain_fraction": (1.0, r"$D_{\mathrm{small}}/D_{\mathrm{dust}}$"),
    }
    if variable in labels:
        return labels[variable]

    # These factors are the same constants used by the analysis PDF renderer.
    if variable == "temperature":
        return 4787.4929974978895, r"$T\;[\mathrm{K}]$"
    if variable == "hydro_w_d":
        return 0.09967741940610811, r"$n\;[\mathrm{cm^{-3}}]$"
    if variable == "hydro_w_e":
        return 477.2049474154028, r"$P/k_B\;[\mathrm{K\,cm^{-3}}]$"
    if variable == "vel_sph_r":
        return 9.777922215131456, r"$v_r\;[\mathrm{km\,s^{-1}}]$"
    if variable == "vel_sph_theta":
        return 9.777922215131456, r"$v_\theta\;[\mathrm{km\,s^{-1}}]$"
    if variable == "vel_sph_phi":
        return 9.777922215131456, r"$v_\phi\;[\mathrm{km\,s^{-1}}]$"
    return 1.0, variable.replace("_", r"\_")


def _axis_quantity_name(variable: str) -> str:
    names = {
        "coord_r": "radius",
        "coord_abscostheta": "costheta",
        "coord_costheta": "signed_costheta",
        "temperature": "temperature",
        "temperature_kelvin": "temperature",
        "hydro_w_d": "density",
        "hydrogen_number_density_cm3": "density",
        "pressure_over_kb_k_cm3": "pressure",
        "cooling_rate_erg_s_cm3": "net_cooling_rate",
        "cooling_time_myr": "net_cooling_time",
        "hydro_w_e": "pressure",
        "vel_sph_r": "radial_velocity",
        "velocity_r_km_s": "radial_velocity",
        "absolute_velocity_r_km_s": "radial_velocity",
        "vel_sph_theta": "polar_velocity",
        "velocity_theta_km_s": "polar_velocity",
        "vel_sph_phi": "azimuthal_velocity",
        "velocity_phi_km_s": "azimuthal_velocity",
        "mach": "mach",
        "radial_mach": "radial_mach",
        "absolute_radial_mach": "radial_mach",
        "sound_speed_km_s": "sound_speed",
        "absolute_z": "absolute_z",
        "cylindrical_radius": "cylindrical_radius",
        "hydro_w_s_00": "metallicity",
        "hydro_w_s_01": "small_dust",
        "hydro_w_s_02": "large_dust",
        "specific_angular_momentum_z": "angular_momentum",
        "specific_angular_momentum_z_kpc_km_s": "angular_momentum",
    }
    return names.get(variable, _safe_component(variable))


def _weight_filename_qualifier(weight: WeightPresentation) -> Optional[str]:
    if weight.filename_qualifier is not None:
        return _safe_component(weight.filename_qualifier)
    if weight.style_key == "volume":
        return None
    if weight.style_key == "mass":
        return "mass_weighted"
    return _safe_component(weight.style_key)


def _legacy_compatibility_weight(
    weight: WeightPresentation, product: str
) -> WeightPresentation:
    if not re.fullmatch(r"output\d+", product):
        return weight
    qualifier = f"legacy_{product}"
    if weight.filename_qualifier:
        qualifier += f"_{weight.filename_qualifier}"
    return replace(weight, filename_qualifier=qualifier)


def _semantic_plot_name(
    components: Sequence[str],
    axes: Sequence[DisplayAxis],
    retained_indices: Iterable[int],
    weight: WeightPresentation,
) -> str:
    retained = set(retained_indices)
    names = [_safe_component(component) for component in components]
    extras = []
    for index, axis in enumerate(axes):
        if index in retained:
            continue
        if axis.raw_variable in {"coord_r", "coord_abscostheta", "coord_costheta"}:
            extras.append(_axis_quantity_name(axis.raw_variable) + "_integrated")
        else:
            extras.append(_axis_quantity_name(axis.raw_variable))
    for extra in extras:
        if extra not in names:
            names.append(extra)
    qualifier = _weight_filename_qualifier(weight)
    if qualifier is not None:
        names.append(qualifier)
    return "_".join(names)


def display_axis(info: Mapping[str, Any]) -> DisplayAxis:
    variable = str(info["variable"])
    factor, label = _axis_factor_and_label(variable)
    return DisplayAxis(
        raw_variable=variable,
        label=label,
        scale=str(info["scale"]),
        edges=np.asarray(info["bin_edges"], dtype=float) * factor,
        centers=np.asarray(info["bin_centers"], dtype=float) * factor,
        linthresh=float(info.get("linthresh", 1.0)) * abs(factor),
    )


def _axis_metadata(axis: DisplayAxis) -> Dict[str, Any]:
    angular_semantics = None
    if axis.raw_variable == "coord_abscostheta":
        angular_semantics = "two hemispheres folded and summed"
    elif axis.raw_variable == "coord_costheta":
        angular_semantics = (
            "north-hemisphere interior; south hemisphere is underflow"
            if axis.edges[0] >= 0.0
            else "signed full-sphere polar cosine"
        )
    return {
        "variable": axis.raw_variable,
        "label": axis.label,
        "scale": axis.scale,
        "display_edges": axis.edges.tolist(),
        "angular_semantics": angular_semantics,
    }


def weight_presentation(header: Mapping[str, Any]) -> WeightPresentation:
    weight = str(header.get("weight", "volume"))
    variable = str(header.get("weight_variable", ""))
    mass_cgs = 3.036951775493658e39
    length_cgs = 3.0856775809623245e21
    time_cgs = 3.15576e15
    mass_msun = mass_cgs / M_SUN_CGS
    mass_rate = mass_cgs / time_cgs / (M_SUN_CGS / YEAR_S)
    energy_rate = mass_cgs * length_cgs**2 / time_cgs**3
    energy = mass_cgs * length_cgs**2 / time_cgs**2

    if weight == "mass":
        return WeightPresentation("mass", r"$M$", r"$M_\odot$", mass_msun, False)
    if weight == "volume":
        return WeightPresentation("volume", r"$V$", r"$\mathrm{kpc^3}$", 1.0, False)
    if variable in {"total_metal_mass", "dust_mass"}:
        label = (
            r"$M_{\mathrm{metal}}$"
            if variable == "total_metal_mass"
            else r"$M_{\mathrm{dust}}$"
        )
        return WeightPresentation(
            "mass",
            label,
            r"$M_\odot$",
            mass_msun,
            False,
            WEIGHT_FILENAME_QUALIFIERS[variable],
        )
    if "mdot" in variable:
        signed = variable in SIGNED_WEIGHT_VARIABLES
        if variable in NEGATIVE_INFLOW_WEIGHT_VARIABLES:
            style = "negative_inflow_mass_flux"
        elif variable.endswith("_in_abs"):
            style = "positive_inflow_mass_flux"
        else:
            style = "signed_net_mass_flux" if signed else "positive_outflow_mass_flux"
        return WeightPresentation(
            style,
            WEIGHT_LABELS.get(variable, r"$\dot{M}$"),
            r"$M_\odot\,\mathrm{kpc\,yr^{-1}}$",
            mass_rate,
            signed,
            WEIGHT_FILENAME_QUALIFIERS.get(variable, _safe_component(variable)),
        )
    if variable == "edot_cool":
        return WeightPresentation(
            "signed_net_energy_flux",
            WEIGHT_LABELS[variable],
            r"$\mathrm{erg\,s^{-1}}$",
            1.0,
            True,
            WEIGHT_FILENAME_QUALIFIERS[variable],
        )
    if "edot" in variable:
        signed = variable in SIGNED_WEIGHT_VARIABLES
        if variable in NEGATIVE_INFLOW_WEIGHT_VARIABLES:
            style = "negative_inflow_energy_flux"
        elif variable.endswith("_in_abs"):
            style = "positive_inflow_energy_flux"
        else:
            style = (
                "signed_net_energy_flux" if signed else "positive_outflow_energy_flux"
            )
        return WeightPresentation(
            style,
            WEIGHT_LABELS.get(variable, r"$\dot{E}$"),
            r"$\mathrm{erg\,kpc\,s^{-1}}$",
            energy_rate,
            signed,
            WEIGHT_FILENAME_QUALIFIERS.get(variable, _safe_component(variable)),
        )
    if "ram" in variable:
        return WeightPresentation(
            "positive_outflow_ram_momentum_flux",
            WEIGHT_LABELS.get(variable, r"$\sum \rho v^2\,dV$"),
            r"$\mathrm{erg}$",
            energy,
            False,
            WEIGHT_FILENAME_QUALIFIERS.get(variable, _safe_component(variable)),
        )
    return WeightPresentation(
        "code",
        r"$W$",
        "code",
        1.0,
        False,
        _safe_component(variable) if variable else "code_weighted",
    )


def _configured_cmap(config: Mapping[str, Any], style_key: str) -> Any:
    entry = config.get("heatmap_weights", {}).get(style_key, {})
    candidates = (
        entry.get("cmap"),
        entry.get("fallback_cmap"),
        CMAP_FALLBACKS.get(style_key, CMAP_FALLBACKS["code"]),
        "Greys" if style_key in NONNEGATIVE_HEATMAP_STYLE_KEYS else None,
    )
    for configured_name in candidates:
        if configured_name is None:
            continue
        configured = _strip_matplotlib_prefix(str(configured_name))
        names = [configured]
        if style_key in NONNEGATIVE_HEATMAP_STYLE_KEYS:
            reversed_name = (
                configured[: -len("_r")]
                if configured.endswith("_r")
                else configured + "_r"
            )
            if reversed_name not in names:
                names.append(reversed_name)
        for name in names:
            try:
                cmap = copy.copy(plt.get_cmap(name))
            except ValueError:
                continue
            if style_key in NONNEGATIVE_HEATMAP_STYLE_KEYS and not np.allclose(
                cmap(0.0)[:3], (1.0, 1.0, 1.0), atol=2.0e-2
            ):
                continue
            policy = config.get("color_policy", {})
            # Empty PDF bins are absence of support, not a physical field value.
            cmap.set_bad(EMPTY_SUPPORT_COLOR)
            if policy.get("under_color") is not None:
                cmap.set_under(policy["under_color"])
            if policy.get("over_color") is not None:
                cmap.set_over(policy["over_color"])
            return cmap
    raise ValueError(f"No available colormap for weight style {style_key!r}")


def _heatmap_cmap(
    config: Mapping[str, Any], style_key: str, norm: colors.Normalize
) -> Any:
    if isinstance(norm, colors.SymLogNorm):
        return _configured_cmap(config, "symlog_diverging")
    return _configured_cmap(config, style_key)


def _weight_norm(
    values: np.ndarray, signed: bool
) -> Tuple[colors.Normalize, np.ma.MaskedArray, Dict[str, Any]]:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        raise ValueError("Marginal contains no finite weights")
    magnitudes = np.abs(finite[finite != 0.0])
    if magnitudes.size == 0:
        raise ValueError("Marginal contains no nonzero weights")
    ordered = np.sort(magnitudes)
    absolute_total = float(np.sum(ordered))
    budget = absolute_total * VISUAL_OMITTED_ABSOLUTE_WEIGHT_FRACTION
    cumulative = np.cumsum(ordered)
    floor_index = int(np.searchsorted(cumulative, budget, side="right"))
    visual_floor = float(ordered[floor_index - 1]) if floor_index else 0.0
    mask = ~np.isfinite(array) | (np.abs(array) <= visual_floor)
    omitted_absolute = float(np.sum(np.abs(array[np.isfinite(array) & mask])))
    plotted = np.ma.masked_where(mask, array)
    visible = plotted.compressed()
    visual_metadata = {
        "visual_floor_absolute": visual_floor,
        "visual_floor_budget_fraction": VISUAL_OMITTED_ABSOLUTE_WEIGHT_FRACTION,
        "omitted_absolute_fraction": (
            omitted_absolute / absolute_total if absolute_total > 0.0 else None
        ),
        "masked_bin_count": int(np.count_nonzero(mask)),
        "visible_bin_count": int(visible.size),
    }
    if signed:
        maximum = float(np.max(np.abs(visible)))
        linthresh = max(maximum * 1.0e-3, np.finfo(float).tiny)
        norm = colors.SymLogNorm(linthresh=linthresh, vmin=-maximum, vmax=maximum)
        return (
            norm,
            plotted,
            {
                "kind": "symlog",
                "vmin": -maximum,
                "vmax": maximum,
                "linthresh": linthresh,
                **visual_metadata,
            },
        )
    positive = visible[visible > 0.0]
    if positive.size == 0:
        raise ValueError("Positive marginal contains no visible positive weights")
    minimum = float(np.min(positive))
    maximum = float(np.max(positive))
    if minimum == maximum:
        maximum = minimum * (1.0 + 1.0e-12)
    norm = colors.LogNorm(vmin=minimum, vmax=maximum)
    return (
        norm,
        plotted,
        {
            "kind": "log",
            "vmin": minimum,
            "vmax": maximum,
            "linthresh": None,
            **visual_metadata,
        },
    )


def _apply_axis_scale(ax: Any, axis: DisplayAxis, coordinate: str) -> None:
    setter = ax.set_xscale if coordinate == "x" else ax.set_yscale
    limiter = ax.set_xlim if coordinate == "x" else ax.set_ylim
    if axis.scale == "log":
        setter("log")
    elif axis.scale == "symlog":
        setter("symlog", linthresh=max(axis.linthresh, np.finfo(float).tiny))
    lower, upper = DISPLAY_AXIS_LIMITS.get(
        axis.raw_variable, (float(axis.edges[0]), float(axis.edges[-1]))
    )
    limiter(float(lower), float(upper))


def _draw_reference_lines(
    ax: Any, axis: DisplayAxis, coordinate: str, config: Mapping[str, Any]
) -> None:
    line = ax.axvline if coordinate == "x" else ax.axhline
    if axis.raw_variable in {"temperature", "temperature_kelvin"}:
        for limits in config.get("temperature_phases_k", {}).values():
            lower = limits[0]
            if lower is not None and axis.edges[0] < float(lower) < axis.edges[-1]:
                line(float(lower), color="white", linewidth=0.7, alpha=0.65)
    elif axis.raw_variable in {"mach", "absolute_radial_mach"}:
        if axis.edges[0] < 1.0 < axis.edges[-1]:
            line(1.0, color="white", linewidth=0.8, alpha=0.75)
    elif axis.raw_variable == "radial_mach":
        for value in (-1.0, 0.0, 1.0):
            if axis.edges[0] < value < axis.edges[-1]:
                line(value, color="white", linewidth=0.8, alpha=0.75)
    elif axis.raw_variable in {
        "coord_costheta",
        "vel_sph_r",
        "vel_sph_theta",
        "vel_sph_phi",
        "velocity_r_km_s",
        "velocity_theta_km_s",
        "velocity_phi_km_s",
        "specific_angular_momentum_z",
        "specific_angular_momentum_z_kpc_km_s",
    }:
        if axis.edges[0] < 0.0 < axis.edges[-1]:
            line(0.0, color="white", linewidth=0.8, alpha=0.75)


def _excluded_summary(data: np.ndarray, interior: np.ndarray) -> Dict[str, Any]:
    absolute_total = float(np.sum(np.abs(data)))
    absolute_interior = float(np.sum(np.abs(interior)))
    return {
        "total_weight_code": float(np.sum(data)),
        "interior_weight_code": float(np.sum(interior)),
        "absolute_total_weight_code": absolute_total,
        "absolute_interior_weight_code": absolute_interior,
        "excluded_absolute_fraction": (
            (absolute_total - absolute_interior) / absolute_total
            if absolute_total > 0.0
            else None
        ),
    }


def _overflow_by_axis(
    data: np.ndarray, axes: Sequence[DisplayAxis]
) -> List[Dict[str, Any]]:
    absolute_total = float(np.sum(np.abs(data)))
    summaries = []
    for index, axis in enumerate(axes):
        under = np.take(data, 0, axis=index)
        over = np.take(data, -1, axis=index)
        under_absolute = float(np.sum(np.abs(under)))
        over_absolute = float(np.sum(np.abs(over)))
        summaries.append(
            {
                "variable": axis.raw_variable,
                "underflow_weight_code": float(np.sum(under)),
                "overflow_weight_code": float(np.sum(over)),
                "underflow_absolute_weight_code": under_absolute,
                "overflow_absolute_weight_code": over_absolute,
                "underflow_absolute_fraction": (
                    under_absolute / absolute_total if absolute_total > 0.0 else None
                ),
                "overflow_absolute_fraction": (
                    over_absolute / absolute_total if absolute_total > 0.0 else None
                ),
            }
        )
    return summaries


def _source_context(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    expected_volume = manifest.get("expected_domain_volume")
    summed_volume = manifest.get("summed_leaf_volume")
    validated_geometry = manifest.get("geometry_to_logical_key_validated") is True
    closure_error = None
    if (
        isinstance(expected_volume, (int, float))
        and isinstance(summed_volume, (int, float))
        and float(expected_volume) > 0.0
    ):
        closure_error = abs(float(summed_volume) - float(expected_volume)) / float(
            expected_volume
        )
    complete = (
        validated_geometry and closure_error is not None and closure_error <= 2.0e-9
    )
    shards = manifest.get("shards_processed")
    if not isinstance(shards, int):
        shards = manifest.get("shards_available")
    annotation = None
    if not complete:
        suffix = f": {shards} shards" if isinstance(shards, int) else ""
        annotation = f"PARTIAL / UNVALIDATED REDUCTION{suffix}"
    return {
        "status": "validated_complete_snapshot"
        if complete
        else "partial_or_unvalidated",
        "annotation": annotation,
        "geometry_to_logical_key_validated": validated_geometry,
        "expected_domain_volume": expected_volume,
        "summed_leaf_volume": summed_volume,
        "relative_volume_closure_error": closure_error,
        "shards_available": manifest.get("shards_available"),
        "shards_processed": manifest.get("shards_processed"),
        "meshblocks_processed": manifest.get("meshblocks_processed"),
        "cells_processed": manifest.get("cells_processed"),
    }


def _decorate_axes(
    ax: Any,
    *,
    source_context: Mapping[str, Any],
) -> None:
    annotation = source_context.get("annotation")
    if annotation:
        ax.text(
            0.98,
            0.98,
            str(annotation),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize="x-small",
            color="#a50f15",
            bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "#a50f15"},
        )


def _save_figure(
    fig: Any,
    base_path: Path,
    formats: Sequence[str],
    *,
    dpi: int,
    bbox_inches: str,
) -> List[str]:
    outputs: List[str] = []
    base_path.parent.mkdir(parents=True, exist_ok=True)
    for extension in formats:
        output = Path(str(base_path) + "." + extension)
        if output.exists():
            raise FileExistsError(
                f"Semantic plot filename collision: {output}. Select fewer "
                "products or make their scientific filename contexts distinct."
            )
        temporary = output.with_name(output.name + f".{os.getpid()}.partial")
        fig.savefig(temporary, format=extension, dpi=dpi, bbox_inches=bbox_inches)
        try:
            os.link(temporary, output)
        except FileExistsError as exc:
            raise FileExistsError(
                f"Semantic plot filename collision: {output}. Select fewer "
                "products or make their scientific filename contexts distinct."
            ) from exc
        finally:
            temporary.unlink(missing_ok=True)
        outputs.append(str(output))
    return outputs


def _marginal_1d(interior: np.ndarray, axis_index: int) -> np.ndarray:
    sum_axes = tuple(index for index in range(interior.ndim) if index != axis_index)
    return np.sum(interior, axis=sum_axes) if sum_axes else np.asarray(interior)


def _marginal_2d(interior: np.ndarray, x_index: int, y_index: int) -> np.ndarray:
    sum_axes = tuple(
        index for index in range(interior.ndim) if index not in {x_index, y_index}
    )
    values = np.sum(interior, axis=sum_axes) if sum_axes else np.asarray(interior)
    retained = [index for index in range(interior.ndim) if index in {x_index, y_index}]
    if retained != [x_index, y_index]:
        values = np.transpose(values)
    return np.asarray(values)


def _angular_bin_selection_weights(
    axis: DisplayAxis, lower_abs_costheta: float, upper_abs_costheta: float
) -> np.ndarray:
    if not 0.0 <= lower_abs_costheta < upper_abs_costheta <= 1.0:
        raise ValueError(
            "Angular sector bounds must satisfy 0 <= lower < upper <= 1"
        )
    edges = np.asarray(axis.edges, dtype=float)
    widths = np.diff(edges)
    if np.any(widths <= 0.0):
        raise ValueError("Angular axis edges must be strictly increasing")
    if axis.raw_variable == "coord_abscostheta":
        overlap = np.maximum(
            0.0,
            np.minimum(edges[1:], upper_abs_costheta)
            - np.maximum(edges[:-1], lower_abs_costheta),
        )
    elif axis.raw_variable == "coord_costheta":
        positive_overlap = np.maximum(
            0.0,
            np.minimum(edges[1:], upper_abs_costheta)
            - np.maximum(edges[:-1], lower_abs_costheta),
        )
        negative_overlap = np.maximum(
            0.0,
            np.minimum(edges[1:], -lower_abs_costheta)
            - np.maximum(edges[:-1], -upper_abs_costheta),
        )
        overlap = positive_overlap + negative_overlap
    else:
        raise ValueError(f"{axis.raw_variable} is not a supported polar-angle axis")
    return overlap / widths


def _bin_range_selection_weights(
    axis: DisplayAxis, lower: Optional[float], upper: Optional[float]
) -> np.ndarray:
    edges = np.asarray(axis.edges, dtype=float)
    widths = np.diff(edges)
    if np.any(widths <= 0.0):
        raise ValueError("Axis edges must be strictly increasing")
    lower_value = float(edges[0]) if lower is None else float(lower)
    upper_value = float(edges[-1]) if upper is None else float(upper)
    overlap = np.maximum(
        0.0,
        np.minimum(edges[1:], upper_value)
        - np.maximum(edges[:-1], lower_value),
    )
    return overlap / widths


def _selected_marginal_2d(
    interior: np.ndarray,
    x_index: int,
    y_index: int,
    selection_index: int,
    selection_weights: np.ndarray,
) -> np.ndarray:
    if selection_index in {x_index, y_index}:
        raise ValueError("Selection axis must not be retained in the marginal")
    shape = [1] * interior.ndim
    shape[selection_index] = selection_weights.size
    selected = interior * selection_weights.reshape(shape)
    return _marginal_2d(selected, x_index, y_index)


def _display_axis_rank(axis: DisplayAxis) -> int:
    variable = axis.raw_variable
    if variable in {"coord_r", "cylindrical_radius"}:
        return 0
    if variable == "absolute_z":
        return 1
    if variable in {"coord_abscostheta", "coord_costheta"}:
        return 3
    if variable in {"temperature", "temperature_kelvin"}:
        return 4
    return 2


def _ordered_pair(
    first_index: int, second_index: int, axes: Sequence[DisplayAxis]
) -> Tuple[int, int]:
    first_key = (_display_axis_rank(axes[first_index]), first_index)
    second_key = (_display_axis_rank(axes[second_index]), second_index)
    return (
        (first_index, second_index)
        if first_key <= second_key
        else (second_index, first_index)
    )


def _redundant_state_geometry_pair(
    x_axis: DisplayAxis, y_axis: DisplayAxis, weight: WeightPresentation
) -> bool:
    variables = {x_axis.raw_variable, y_axis.raw_variable}
    return weight.style_key in {"mass", "volume"} and variables in (
        {"coord_r", "coord_abscostheta"},
        {"coord_r", "coord_costheta"},
    )


def _angular_axis_index(axes: Sequence[DisplayAxis]) -> Optional[int]:
    matches = [
        index
        for index, axis in enumerate(axes)
        if axis.raw_variable in {"coord_abscostheta", "coord_costheta"}
    ]
    if len(matches) > 1:
        raise ValueError("A PDF product may contain at most one polar-angle axis")
    return matches[0] if matches else None


def _radial_axis_index(axes: Sequence[DisplayAxis]) -> Optional[int]:
    matches = [
        index for index, axis in enumerate(axes) if axis.raw_variable == "coord_r"
    ]
    if len(matches) > 1:
        raise ValueError("A PDF product may contain at most one radial axis")
    return matches[0] if matches else None


def _sectorized_radial_pair(
    x_index: int,
    y_index: int,
    *,
    radial_index: Optional[int],
    angular_index: Optional[int],
) -> bool:
    return (
        radial_index is not None
        and angular_index is not None
        and radial_index in {x_index, y_index}
        and angular_index not in {x_index, y_index}
    )


def _replaced_by_sectorized_radial_pair(
    x_index: int,
    y_index: int,
    *,
    radial_index: Optional[int],
    angular_index: Optional[int],
) -> bool:
    return (
        radial_index is not None
        and angular_index is not None
        and angular_index in {x_index, y_index}
        and radial_index not in {x_index, y_index}
    )


def _is_radial_transport(weight_variable: str) -> bool:
    return weight_variable.startswith(("mdot_sph", "edot_sph")) or weight_variable == (
        "radial_ram_out"
    )


def _crossing_presentation(
    weight: WeightPresentation, weight_variable: str
) -> WeightPresentation:
    # The shell-integrated transport products are displayed as physical rates,
    # not as their underlying volume-integral bookkeeping expression.
    label = weight.label
    if weight_variable.startswith("mdot_sph"):
        unit = r"$M_\odot\,\mathrm{yr^{-1}}$"
    elif weight_variable.startswith("edot_sph"):
        unit = r"$\mathrm{erg\,s^{-1}}$"
    elif weight_variable == "radial_ram_out":
        label = r"$\dot{p}_{\mathrm{out}}$"
        unit = r"$\mathrm{erg\,kpc^{-1}}$"
    else:
        unit = r"$\mathrm{erg\,kpc^{-1}}$"
    return WeightPresentation(
        weight.style_key,
        label,
        unit,
        weight.factor,
        weight.signed,
        weight.filename_qualifier,
    )


def _shell_volume_fraction_presentation() -> WeightPresentation:
    return WeightPresentation(
        "volume_fraction", r"$V_{\mathrm{bin}}/V_{\mathrm{shell}}$", "", 1.0, False
    )


def _wedge_volume_fraction_presentation() -> WeightPresentation:
    return WeightPresentation(
        "volume_fraction", r"$V_{\mathrm{bin}}/V_{\mathrm{wedge}}$", "", 1.0, False
    )


def _weight_axis_label(weight: WeightPresentation) -> str:
    return f"{weight.label} [{weight.unit}]" if weight.unit else weight.label


def _mathtext_body(label: str) -> str:
    if label.startswith("$") and label.endswith("$"):
        return label[1:-1]
    return label


def _differential_axis_symbol(axis: DisplayAxis) -> str:
    symbol = _mathtext_body(axis.label)
    if r"\;[" in symbol:
        symbol = symbol.split(r"\;[", 1)[0]
    return symbol


def _differential_bin_widths(axis: DisplayAxis) -> np.ndarray:
    edges = np.asarray(axis.edges, dtype=float)
    if axis.scale == "log":
        if np.any(edges <= 0.0):
            raise ValueError(
                f"Log-scaled {axis.raw_variable} axis must have positive edges"
            )
        return np.diff(np.log10(edges))
    return np.diff(edges)


def _differential_profile_axis_label(
    weight: WeightPresentation, axis: DisplayAxis
) -> str:
    numerator = _mathtext_body(weight.label)
    denominator = _differential_axis_symbol(axis)
    differential = r"\mathrm{d}"
    if axis.scale == "log":
        return (
            rf"$\frac{{{differential}{numerator}}}"
            rf"{{{differential}\log_{{10}} {denominator}}}$"
        )
    return rf"$\frac{{{differential}{numerator}}}{{{differential}{denominator}}}$"


def _radial_bin_width_normalize(
    values: np.ndarray,
    axes: Sequence[DisplayAxis],
    *,
    radial_axis_index: int,
) -> np.ndarray:
    widths = np.diff(axes[radial_axis_index].edges)
    shape = [1] * values.ndim
    shape[radial_axis_index] = widths.size
    return values / widths.reshape(shape)


def _spherical_shell_volume_normalize(
    values: np.ndarray,
    axes: Sequence[DisplayAxis],
    *,
    radial_axis_index: int,
    shell_fraction: float = 1.0,
) -> np.ndarray:
    if not math.isfinite(shell_fraction) or not 0.0 < shell_fraction <= 1.0:
        raise ValueError("Spherical-shell fraction must be within (0, 1]")
    edges = np.asarray(axes[radial_axis_index].edges, dtype=float)
    shell_volumes = (
        shell_fraction
        * (4.0 * math.pi / 3.0)
        * (edges[1:] ** 3 - edges[:-1] ** 3)
    )
    if np.any(~np.isfinite(shell_volumes)) or np.any(shell_volumes <= 0.0):
        raise ValueError("Radial edges do not define positive spherical-shell volumes")
    shape = [1] * values.ndim
    shape[radial_axis_index] = shell_volumes.size
    return values / shell_volumes.reshape(shape)


def _radial_angular_wedge_volume_normalize(
    values: np.ndarray,
    angular_axis: DisplayAxis,
    radial_axis: DisplayAxis,
    *,
    angular_axis_index: int,
    radial_bin_index: int,
) -> np.ndarray:
    radial_edges = np.asarray(radial_axis.edges, dtype=float)
    shell_volume = (4.0 * math.pi / 3.0) * (
        radial_edges[radial_bin_index + 1] ** 3 - radial_edges[radial_bin_index] ** 3
    )
    angular_widths = np.diff(np.asarray(angular_axis.edges, dtype=float))
    if angular_axis.raw_variable == "coord_abscostheta":
        wedge_fractions = angular_widths
    elif angular_axis.raw_variable == "coord_costheta":
        wedge_fractions = 0.5 * angular_widths
    else:
        raise ValueError(f"{angular_axis.raw_variable} is not a polar-angle axis")
    wedge_volumes = shell_volume * wedge_fractions
    if np.any(~np.isfinite(wedge_volumes)) or np.any(wedge_volumes <= 0.0):
        raise ValueError("Radial/angular bins do not define positive wedge volumes")
    shape = [1] * values.ndim
    shape[angular_axis_index] = wedge_volumes.size
    return values / wedge_volumes.reshape(shape)


def _prepare_heatmap_values(
    values: np.ndarray,
    x_axis: DisplayAxis,
    y_axis: DisplayAxis,
    weight: WeightPresentation,
    weight_variable: str,
    *,
    shell_fraction: float = 1.0,
) -> Tuple[np.ndarray, WeightPresentation, Optional[str], Optional[str]]:
    panel_weight = weight
    transport_normalization = None
    radial_normalization = None
    radial_axis_index = None
    if x_axis.raw_variable == "coord_r":
        radial_axis_index = 0
    elif y_axis.raw_variable == "coord_r":
        radial_axis_index = 1
    if radial_axis_index is not None:
        if weight.style_key == "volume":
            values = _spherical_shell_volume_normalize(
                values,
                [x_axis, y_axis],
                radial_axis_index=radial_axis_index,
                shell_fraction=shell_fraction,
            )
            panel_weight = _shell_volume_fraction_presentation()
            radial_normalization = "divided_by_angular_sector_shell_volume"
        elif _is_radial_transport(weight_variable):
            values = _radial_bin_width_normalize(
                values, [x_axis, y_axis], radial_axis_index=radial_axis_index
            )
            panel_weight = _crossing_presentation(weight, weight_variable)
            transport_normalization = "divided_by_radial_bin_width"
    return (
        np.asarray(values, dtype=float) * panel_weight.factor,
        panel_weight,
        transport_normalization,
        radial_normalization,
    )


def _angular_sector_bounds(
    config: Mapping[str, Any], sector: str
) -> Tuple[float, float]:
    values = config["angular_selections"][sector]["abs_costheta_range"]
    return float(values[0]), float(values[1])


def _weighted_mean_median(
    weights: np.ndarray, physical_centers: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    values = np.asarray(weights, dtype=float)
    centers = np.asarray(physical_centers, dtype=float)
    if values.ndim != 2 or values.shape[1] != centers.size:
        raise ValueError("Profile weights must have shape (profile_bin, physical_bin)")
    if np.any(values < 0.0):
        raise ValueError("Conditional means and medians require nonnegative weights")
    means = np.full(values.shape[0], np.nan)
    medians = np.full(values.shape[0], np.nan)
    for index, row in enumerate(values):
        valid = np.isfinite(row) & np.isfinite(centers) & (row > 0.0)
        if not np.any(valid):
            continue
        valid_centers = centers[valid]
        valid_weights = row[valid]
        order = np.argsort(valid_centers)
        valid_centers = valid_centers[order]
        valid_weights = valid_weights[order]
        total = float(np.sum(valid_weights))
        means[index] = float(np.sum(valid_weights * valid_centers) / total)
        cumulative = np.cumsum(valid_weights)
        medians[index] = float(
            np.interp(0.5 * total, cumulative, valid_centers)
        )
    return means, medians


def _finite_list(values: np.ndarray) -> List[Optional[float]]:
    return [
        float(value) if np.isfinite(value) else None
        for value in np.asarray(values, dtype=float)
    ]


def _plot_conditional_profile_summary(
    interior: np.ndarray,
    axes: Sequence[DisplayAxis],
    radial_index: int,
    angular_index: int,
    physical_index: int,
    *,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    radial_axis = axes[radial_index]
    angular_axis = axes[angular_index]
    physical_axis = axes[physical_index]
    fig, subplot_axes = plt.subplots(
        2, 1, figsize=(9.0, 10.5), squeeze=False, constrained_layout=True
    )
    radial_ax, angular_ax = subplot_axes[:, 0]
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    fig.suptitle(time_title, fontfamily="monospace")
    _, _, apply_axes_style = _analysis_api(str(analysis_root))

    sector_colors = dict(
        zip(ANGULAR_SECTOR_ORDER, plt.get_cmap("tab10").colors, strict=False)
    )
    radial_series = []
    for sector in ANGULAR_SECTOR_ORDER:
        lower, upper = _angular_sector_bounds(config, sector)
        selection_weights = _angular_bin_selection_weights(
            angular_axis, lower, upper
        )
        marginal = _selected_marginal_2d(
            interior,
            radial_index,
            physical_index,
            angular_index,
            selection_weights,
        )
        mean, median = _weighted_mean_median(marginal, physical_axis.centers)
        color = sector_colors[sector]
        radial_ax.plot(
            radial_axis.centers,
            mean,
            color=color,
            linewidth=1.7,
            label=ANGULAR_SECTOR_TITLES[sector],
        )
        radial_ax.plot(
            radial_axis.centers,
            median,
            color=color,
            linestyle="--",
            linewidth=1.5,
        )
        radial_series.append(
            {
                "sector": sector,
                "abs_costheta_range": [lower, upper],
                "mean": _finite_list(mean),
                "median": _finite_list(median),
            }
        )
    radial_ax.set_xlabel(radial_axis.label)
    radial_ax.set_ylabel(physical_axis.label)
    radial_ax.set_title("Angular selections", fontsize="small")
    _apply_axis_scale(radial_ax, radial_axis, "x")
    _apply_axis_scale(radial_ax, physical_axis, "y")
    _decorate_axes(radial_ax, source_context=source_context)
    apply_axes_style(radial_ax, config["axes_style"])
    sector_legend = radial_ax.legend(loc="best", fontsize="small", ncol=2)
    radial_ax.add_artist(sector_legend)
    radial_ax.legend(
        handles=[
            Line2D([0], [0], color="#252525", linewidth=1.7, label="mean"),
            Line2D(
                [0],
                [0],
                color="#252525",
                linestyle="--",
                linewidth=1.5,
                label="median",
            ),
        ],
        loc="lower left",
        fontsize="small",
    )

    shell_cmap = plt.get_cmap("viridis")
    angular_series = []
    for index, target_radius in enumerate(FIXED_RADIAL_SHELL_TARGETS_KPC):
        selected_index = int(
            np.argmin(np.abs(np.log(radial_axis.centers / target_radius)))
        )
        selection_weights = np.zeros(radial_axis.centers.size, dtype=float)
        selection_weights[selected_index] = 1.0
        marginal = _selected_marginal_2d(
            interior,
            angular_index,
            physical_index,
            radial_index,
            selection_weights,
        )
        mean, median = _weighted_mean_median(marginal, physical_axis.centers)
        color = shell_cmap(index / (len(FIXED_RADIAL_SHELL_TARGETS_KPC) - 1))
        angular_ax.plot(
            angular_axis.centers,
            mean,
            color=color,
            linewidth=1.7,
            label=rf"$r \simeq {target_radius:g}\;\mathrm{{kpc}}$",
        )
        angular_ax.plot(
            angular_axis.centers,
            median,
            color=color,
            linestyle="--",
            linewidth=1.5,
        )
        angular_series.append(
            {
                "target_radius_kpc": target_radius,
                "selected_radial_bin_index": selected_index,
                "selected_radial_center_kpc": float(
                    radial_axis.centers[selected_index]
                ),
                "selected_radial_edges_kpc": [
                    float(radial_axis.edges[selected_index]),
                    float(radial_axis.edges[selected_index + 1]),
                ],
                "mean": _finite_list(mean),
                "median": _finite_list(median),
            }
        )
    angular_ax.set_xlabel(angular_axis.label)
    angular_ax.set_ylabel(physical_axis.label)
    angular_ax.set_title("Fixed radial shells", fontsize="small")
    _apply_axis_scale(angular_ax, angular_axis, "x")
    _apply_axis_scale(angular_ax, physical_axis, "y")
    _decorate_axes(angular_ax, source_context=source_context)
    apply_axes_style(angular_ax, config["axes_style"])
    angular_ax.legend(loc="best", fontsize="small", ncol=2)

    paths = _save_figure(fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
    return {
        "kind": "conditional_mean_median_profiles",
        "paths": paths,
        "radial_axis": radial_axis.raw_variable,
        "angular_axis": angular_axis.raw_variable,
        "physical_axis": physical_axis.raw_variable,
        "statistics": ["mean", "median"],
        "radial_series": radial_series,
        "angular_series": angular_series,
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _plot_heatmap(
    values: np.ndarray,
    x_axis: DisplayAxis,
    y_axis: DisplayAxis,
    weight: WeightPresentation,
    weight_variable: str,
    excluded: Mapping[str, Any],
    *,
    product: str,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    _, add_matched_colorbar, apply_axes_style = _analysis_api(str(analysis_root))
    figure_size = config["figure_sizes_inches"]["pdf_heatmap"]
    fig, ax = plt.subplots(figsize=figure_size)
    scaled, panel_weight, transport_normalization, radial_normalization = (
        _prepare_heatmap_values(values, x_axis, y_axis, weight, weight_variable)
    )
    ax.set_xlabel(x_axis.label)
    ax.set_ylabel(y_axis.label)
    _apply_axis_scale(ax, x_axis, "x")
    _apply_axis_scale(ax, y_axis, "y")
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    ax.set_title(time_title, fontfamily="monospace")
    _decorate_axes(ax, source_context=source_context)
    apply_axes_style(ax, config["axes_style"])
    status = "ok"
    norm_metadata = None
    cmap = None
    if np.any(scaled != 0.0):
        norm, plotted, norm_metadata = _weight_norm(scaled, panel_weight.signed)
        cmap = _heatmap_cmap(config, panel_weight.style_key, norm)
        artist = ax.pcolormesh(
            x_axis.edges,
            y_axis.edges,
            plotted.T,
            shading="flat",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        _draw_reference_lines(ax, x_axis, "x", config)
        _draw_reference_lines(ax, y_axis, "y", config)
        if {
            x_axis.raw_variable,
            y_axis.raw_variable,
        } == {"absolute_velocity_r_km_s", "sound_speed_km_s"}:
            lower = max(float(x_axis.edges[0]), float(y_axis.edges[0]))
            upper = min(float(x_axis.edges[-1]), float(y_axis.edges[-1]))
            ax.plot(
                [lower, upper],
                [lower, upper],
                color="#252525",
                linestyle=":",
                linewidth=1.2,
                label=r"$|v_r|=c_s$",
            )
            positive_weights = np.abs(np.asarray(scaled, dtype=float))
            _, median_y = _weighted_mean_median(
                positive_weights, y_axis.centers
            )
            _, median_x = _weighted_mean_median(
                positive_weights.T, x_axis.centers
            )
            ax.plot(
                x_axis.centers,
                median_y,
                color="#111111",
                linewidth=1.3,
                label=rf"median {y_axis.label}",
            )
            ax.plot(
                median_x,
                y_axis.centers,
                color="#555555",
                linestyle="--",
                linewidth=1.3,
                label=rf"median {x_axis.label}",
            )
            ax.legend(loc="best", fontsize="x-small")
        add_matched_colorbar(
            fig,
            ax,
            artist,
            label=_weight_axis_label(panel_weight),
            style=config["axes_style"],
        )
    else:
        status = "no_interior_support"
        ax.set_facecolor(EMPTY_SUPPORT_COLOR)
        excluded_fraction = excluded.get("excluded_absolute_fraction")
        detail = (
            "100% OF WEIGHT IS IN UNDER/OVERFLOW"
            if excluded_fraction is not None
            and float(excluded_fraction) >= 1.0 - 1.0e-12
            else "ZERO NET WEIGHT IN THIS MARGINAL"
        )
        ax.text(
            0.5,
            0.5,
            "NO JOINT INTERIOR SUPPORT\n" + detail,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize="small",
            fontweight="bold",
        )
    paths = _save_figure(fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
    return {
        "kind": "two_dimensional_marginal",
        "paths": paths,
        "x_axis": x_axis.raw_variable,
        "y_axis": y_axis.raw_variable,
        "status": status,
        "minimum": float(np.nanmin(scaled)),
        "maximum": float(np.nanmax(scaled)),
        "sum": float(np.nansum(scaled)),
        "transport_normalization": transport_normalization,
        "radial_normalization": radial_normalization,
        "weight_label": panel_weight.label,
        "weight_unit": panel_weight.unit,
        "cmap": cmap.name if cmap is not None else None,
        "norm": norm_metadata,
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _style_profile_y_axis(ax: Any, values: Sequence[np.ndarray], signed: bool) -> None:
    finite = np.concatenate(
        [np.asarray(value, dtype=float)[np.isfinite(value)] for value in values]
    )
    nonzero = finite[finite != 0.0]
    if nonzero.size == 0:
        return
    if signed or np.any(nonzero < 0.0):
        maximum = float(np.max(np.abs(nonzero)))
        ax.set_yscale(
            "symlog", linthresh=max(maximum * 1.0e-4, np.finfo(float).tiny)
        )
    else:
        ax.set_yscale("log")


def _plot_selected_profiles_2d(
    interior: np.ndarray,
    axes: Sequence[DisplayAxis],
    profile_index: int,
    selection_index: int,
    selections: Sequence[Tuple[str, str, np.ndarray]],
    weight: WeightPresentation,
    weight_variable: str,
    *,
    kind: str,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    _, _, apply_axes_style = _analysis_api(str(analysis_root))
    profile_axis = axes[profile_index]
    selection_axis = axes[selection_index]
    values = _marginal_2d(interior, profile_index, selection_index)
    panel_weight = weight
    transport_normalization = None
    series = []
    plotted_values = []
    for key, label, selection_weights in selections:
        profile = np.sum(
            values * np.asarray(selection_weights, dtype=float)[None, :],
            axis=1,
        )
        if (
            profile_axis.raw_variable == "coord_r"
            and _is_radial_transport(weight_variable)
        ):
            profile = profile / np.diff(profile_axis.edges)
            panel_weight = _crossing_presentation(weight, weight_variable)
            transport_normalization = "divided_by_radial_bin_width"
        scaled = np.asarray(profile, dtype=float) * panel_weight.factor
        plotted_values.append(scaled)
        series.append(
            {
                "selection": key,
                "label": label,
                "values": _finite_list(scaled),
            }
        )

    fig, ax = plt.subplots(figsize=(9.5, 6.5), constrained_layout=True)
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    ax.set_title(time_title, fontfamily="monospace")
    colors_cycle = plt.get_cmap("tab10").colors
    for index, ((_, label, _), profile) in enumerate(
        zip(selections, plotted_values, strict=True)
    ):
        visible = np.where(profile != 0.0, profile, np.nan)
        ax.plot(
            profile_axis.centers,
            visible,
            linewidth=1.7,
            color=colors_cycle[index % len(colors_cycle)],
            label=label,
        )
    ax.set_xlabel(profile_axis.label)
    ax.set_ylabel(_weight_axis_label(panel_weight))
    _apply_axis_scale(ax, profile_axis, "x")
    _style_profile_y_axis(ax, plotted_values, panel_weight.signed)
    _decorate_axes(ax, source_context=source_context)
    apply_axes_style(ax, config["axes_style"])
    ax.legend(loc="best", fontsize="small")
    paths = _save_figure(
        fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches
    )
    plt.close(fig)
    return {
        "kind": kind,
        "paths": paths,
        "profile_axis": profile_axis.raw_variable,
        "selection_axis": selection_axis.raw_variable,
        "series": series,
        "transport_normalization": transport_normalization,
        "weight_label": panel_weight.label,
        "weight_unit": panel_weight.unit,
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _plot_fixed_primary_profiles_2d(
    interior: np.ndarray,
    axes: Sequence[DisplayAxis],
    primary_index: int,
    distribution_index: int,
    weight: WeightPresentation,
    weight_variable: str,
    *,
    kind: str,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    _, _, apply_axes_style = _analysis_api(str(analysis_root))
    primary_axis = axes[primary_index]
    distribution_axis = axes[distribution_index]
    values = _marginal_2d(interior, primary_index, distribution_index)
    panel_weight = weight
    transport_normalization = None
    distribution_normalization = (
        "divided_by_log10_distribution_bin_width"
        if distribution_axis.scale == "log"
        else "divided_by_distribution_bin_width"
    )
    distribution_widths = _differential_bin_widths(distribution_axis)
    series = []
    plotted_values = []
    for target in FIXED_RADIAL_SHELL_TARGETS_KPC:
        selected_index = int(
            np.argmin(np.abs(np.log(primary_axis.centers / target)))
        )
        distribution = (
            np.asarray(values[selected_index, :], dtype=float) / distribution_widths
        )
        if (
            primary_axis.raw_variable == "coord_r"
            and _is_radial_transport(weight_variable)
        ):
            distribution = distribution / np.diff(primary_axis.edges)[selected_index]
            panel_weight = _crossing_presentation(weight, weight_variable)
            transport_normalization = "divided_by_radial_bin_width"
        scaled = distribution * panel_weight.factor
        plotted_values.append(scaled)
        series.append(
            {
                "target_radius_kpc": target,
                "selected_primary_bin_index": selected_index,
                "selected_primary_center": float(
                    primary_axis.centers[selected_index]
                ),
                "values": _finite_list(scaled),
            }
        )

    fig, ax = plt.subplots(figsize=(9.5, 6.5), constrained_layout=True)
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    ax.set_title(time_title, fontfamily="monospace")
    cmap = plt.get_cmap("viridis")
    for index, (target, distribution) in enumerate(
        zip(FIXED_RADIAL_SHELL_TARGETS_KPC, plotted_values, strict=True)
    ):
        visible = np.where(distribution != 0.0, distribution, np.nan)
        ax.plot(
            distribution_axis.centers,
            visible,
            linewidth=1.5,
            color=cmap(index / (len(FIXED_RADIAL_SHELL_TARGETS_KPC) - 1)),
            label=rf"$r \simeq {target:g}\;\mathrm{{kpc}}$",
        )
    ax.set_xlabel(distribution_axis.label)
    ax.set_ylabel(_differential_profile_axis_label(panel_weight, distribution_axis))
    _apply_axis_scale(ax, distribution_axis, "x")
    _draw_reference_lines(ax, distribution_axis, "x", config)
    _style_profile_y_axis(ax, plotted_values, panel_weight.signed)
    _decorate_axes(ax, source_context=source_context)
    apply_axes_style(ax, config["axes_style"])
    ax.legend(loc="best", fontsize="small", ncol=2)
    paths = _save_figure(
        fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches
    )
    plt.close(fig)
    return {
        "kind": kind,
        "paths": paths,
        "primary_axis": primary_axis.raw_variable,
        "distribution_axis": distribution_axis.raw_variable,
        "series": series,
        "transport_normalization": transport_normalization,
        "distribution_normalization": distribution_normalization,
        "weight_label": panel_weight.label,
        "weight_unit": panel_weight.unit,
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _plot_integrated_profile_1d(
    interior: np.ndarray,
    axes: Sequence[DisplayAxis],
    profile_index: int,
    weight: WeightPresentation,
    *,
    kind: str,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    _, _, apply_axes_style = _analysis_api(str(analysis_root))
    profile_axis = axes[profile_index]
    values = _marginal_1d(interior, profile_index) * weight.factor
    fig, ax = plt.subplots(figsize=(9.5, 6.5), constrained_layout=True)
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    ax.set_title(time_title, fontfamily="monospace")
    ax.plot(
        profile_axis.centers,
        np.where(values != 0.0, values, np.nan),
        color="#252525",
        linewidth=1.7,
    )
    ax.set_xlabel(profile_axis.label)
    ax.set_ylabel(_weight_axis_label(weight))
    _apply_axis_scale(ax, profile_axis, "x")
    _style_profile_y_axis(ax, [values], weight.signed)
    _decorate_axes(ax, source_context=source_context)
    apply_axes_style(ax, config["axes_style"])
    paths = _save_figure(
        fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches
    )
    plt.close(fig)
    return {
        "kind": kind,
        "paths": paths,
        "profile_axis": profile_axis.raw_variable,
        "values": _finite_list(values),
        "weight_label": weight.label,
        "weight_unit": weight.unit,
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _plot_angular_sector_heatmaps(
    interior: np.ndarray,
    axes: Sequence[DisplayAxis],
    x_index: int,
    y_index: int,
    angular_index: int,
    weight: WeightPresentation,
    weight_variable: str,
    *,
    product: str,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    _, _, apply_axes_style = _analysis_api(str(analysis_root))
    x_axis = axes[x_index]
    y_axis = axes[y_index]
    angular_axis = axes[angular_index]
    scaled_panels = []
    panel_metadata = []
    panel_weight: Optional[WeightPresentation] = None
    transport_normalization = None
    radial_normalization = None
    for sector in ANGULAR_SECTOR_ORDER:
        lower, upper = _angular_sector_bounds(config, sector)
        selection_weights = _angular_bin_selection_weights(angular_axis, lower, upper)
        marginal = _selected_marginal_2d(
            interior, x_index, y_index, angular_index, selection_weights
        )
        scaled, current_weight, current_transport, current_radial = (
            _prepare_heatmap_values(
                marginal,
                x_axis,
                y_axis,
                weight,
                weight_variable,
                shell_fraction=upper - lower,
            )
        )
        if panel_weight is not None and current_weight != panel_weight:
            raise ValueError("Angular-sector panels do not share weight presentation")
        panel_weight = current_weight
        transport_normalization = current_transport
        radial_normalization = current_radial
        scaled_panels.append(scaled)
        panel_metadata.append(
            {
                "sector": sector,
                "title": ANGULAR_SECTOR_TITLES[sector],
                "abs_costheta_range": [lower, upper],
                "shell_fraction": upper - lower,
                "fractional_boundary_bins": bool(
                    np.any((selection_weights > 0.0) & (selection_weights < 1.0))
                ),
                "minimum": float(np.nanmin(scaled)),
                "maximum": float(np.nanmax(scaled)),
                "sum": float(np.nansum(scaled)),
            }
        )
    assert panel_weight is not None
    stacked = np.stack(scaled_panels)
    status = "ok"
    norm_metadata = None
    plotted_panels: Sequence[np.ma.MaskedArray]
    norm: Optional[colors.Normalize] = None
    if np.any(stacked != 0.0):
        norm, plotted, norm_metadata = _weight_norm(stacked, panel_weight.signed)
        plotted_panels = [plotted[index] for index in range(len(ANGULAR_SECTOR_ORDER))]
    else:
        status = "no_interior_support"
        plotted_panels = [np.ma.masked_all_like(values) for values in scaled_panels]

    fig, subplot_axes = plt.subplots(
        2,
        2,
        figsize=(12.5, 9.5),
        squeeze=False,
        constrained_layout=True,
    )
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    fig.suptitle(time_title, fontfamily="monospace")
    artist = None
    cmap = (
        _heatmap_cmap(config, panel_weight.style_key, norm)
        if norm is not None
        else None
    )
    for index, sector in enumerate(ANGULAR_SECTOR_ORDER):
        ax = subplot_axes.flat[index]
        ax.set_title(ANGULAR_SECTOR_TITLES[sector], fontsize="small")
        ax.set_xlabel(x_axis.label)
        ax.set_ylabel(y_axis.label)
        _apply_axis_scale(ax, x_axis, "x")
        _apply_axis_scale(ax, y_axis, "y")
        _decorate_axes(ax, source_context=source_context)
        apply_axes_style(ax, config["axes_style"])
        if norm is None:
            ax.set_facecolor(EMPTY_SUPPORT_COLOR)
            ax.text(
                0.5,
                0.5,
                "NO JOINT INTERIOR SUPPORT",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize="small",
                fontweight="bold",
            )
            continue
        artist = ax.pcolormesh(
            x_axis.edges,
            y_axis.edges,
            plotted_panels[index].T,
            shading="flat",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        _draw_reference_lines(ax, x_axis, "x", config)
        _draw_reference_lines(ax, y_axis, "y", config)
    if artist is not None:
        colorbar = fig.colorbar(
            artist,
            ax=list(subplot_axes.flat),
            fraction=0.025,
            pad=0.02,
        )
        colorbar.set_label(_weight_axis_label(panel_weight))
        colorbar.ax.tick_params(
            which="both", direction=config["axes_style"]["tick_direction"]
        )
    paths = _save_figure(fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
    return {
        "kind": "angular_sector_radial_marginals",
        "paths": paths,
        "x_axis": x_axis.raw_variable,
        "y_axis": y_axis.raw_variable,
        "selection_axis": angular_axis.raw_variable,
        "status": status,
        "panels": panel_metadata,
        "transport_normalization": transport_normalization,
        "radial_normalization": radial_normalization,
        "weight_label": panel_weight.label,
        "weight_unit": panel_weight.unit,
        "cmap": cmap.name if cmap is not None else None,
        "norm": norm_metadata,
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _plot_fixed_radial_shell_heatmaps(
    interior: np.ndarray,
    axes: Sequence[DisplayAxis],
    angular_index: int,
    physical_index: int,
    radial_index: int,
    weight: WeightPresentation,
    weight_variable: str,
    *,
    product: str,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    _, _, apply_axes_style = _analysis_api(str(analysis_root))
    angular_axis = axes[angular_index]
    physical_axis = axes[physical_index]
    radial_axis = axes[radial_index]
    scaled_panels = []
    panel_metadata = []
    panel_weight = weight
    transport_normalization = None
    radial_normalization = None
    for target_radius in FIXED_RADIAL_SHELL_TARGETS_KPC:
        distances = np.abs(np.log(radial_axis.centers / target_radius))
        selected_index = int(np.argmin(distances))
        selection_weights = np.zeros(radial_axis.centers.size, dtype=float)
        selection_weights[selected_index] = 1.0
        marginal = _selected_marginal_2d(
            interior,
            angular_index,
            physical_index,
            radial_index,
            selection_weights,
        )
        if weight.style_key == "volume":
            marginal = _radial_angular_wedge_volume_normalize(
                marginal,
                angular_axis,
                radial_axis,
                angular_axis_index=0,
                radial_bin_index=selected_index,
            )
            panel_weight = _wedge_volume_fraction_presentation()
            radial_normalization = "divided_by_radial_angular_wedge_volume"
        elif _is_radial_transport(weight_variable):
            radial_width = np.diff(radial_axis.edges)[selected_index]
            marginal = marginal / radial_width
            panel_weight = _crossing_presentation(weight, weight_variable)
            transport_normalization = "divided_by_radial_bin_width"
        scaled = np.asarray(marginal, dtype=float) * panel_weight.factor
        scaled_panels.append(scaled)
        panel_metadata.append(
            {
                "target_radius_kpc": target_radius,
                "selected_radial_bin_index": selected_index,
                "selected_radial_center_kpc": float(radial_axis.centers[selected_index]),
                "selected_radial_edges_kpc": [
                    float(radial_axis.edges[selected_index]),
                    float(radial_axis.edges[selected_index + 1]),
                ],
                "minimum": float(np.nanmin(scaled)),
                "maximum": float(np.nanmax(scaled)),
                "sum": float(np.nansum(scaled)),
            }
        )

    stacked = np.stack(scaled_panels)
    status = "ok"
    norm_metadata = None
    plotted_panels: Sequence[np.ma.MaskedArray]
    norm: Optional[colors.Normalize] = None
    if np.any(stacked != 0.0):
        norm, plotted, norm_metadata = _weight_norm(stacked, panel_weight.signed)
        plotted_panels = [
            plotted[index] for index in range(len(FIXED_RADIAL_SHELL_TARGETS_KPC))
        ]
    else:
        status = "no_interior_support"
        plotted_panels = [np.ma.masked_all_like(values) for values in scaled_panels]

    fig, subplot_axes = plt.subplots(
        3,
        2,
        figsize=(12.5, 13.0),
        squeeze=False,
        constrained_layout=True,
    )
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    fig.suptitle(time_title, fontfamily="monospace")
    artist = None
    cmap = (
        _heatmap_cmap(config, panel_weight.style_key, norm)
        if norm is not None
        else None
    )
    for index, target_radius in enumerate(FIXED_RADIAL_SHELL_TARGETS_KPC):
        ax = subplot_axes.flat[index]
        ax.set_title(rf"$r \simeq {target_radius:g}\;\mathrm{{kpc}}$", fontsize="small")
        ax.set_xlabel(angular_axis.label)
        ax.set_ylabel(physical_axis.label)
        _apply_axis_scale(ax, angular_axis, "x")
        _apply_axis_scale(ax, physical_axis, "y")
        _decorate_axes(ax, source_context=source_context)
        apply_axes_style(ax, config["axes_style"])
        if norm is None:
            ax.set_facecolor(EMPTY_SUPPORT_COLOR)
            ax.text(
                0.5,
                0.5,
                "NO JOINT INTERIOR SUPPORT",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize="small",
                fontweight="bold",
            )
            continue
        artist = ax.pcolormesh(
            angular_axis.edges,
            physical_axis.edges,
            plotted_panels[index].T,
            shading="flat",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        _draw_reference_lines(ax, angular_axis, "x", config)
        _draw_reference_lines(ax, physical_axis, "y", config)
    if artist is not None:
        colorbar = fig.colorbar(
            artist,
            ax=list(subplot_axes.flat),
            fraction=0.025,
            pad=0.02,
        )
        colorbar.set_label(_weight_axis_label(panel_weight))
        colorbar.ax.tick_params(
            which="both", direction=config["axes_style"]["tick_direction"]
        )
    paths = _save_figure(fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
    return {
        "kind": "fixed_radial_shell_angular_marginals",
        "paths": paths,
        "x_axis": angular_axis.raw_variable,
        "y_axis": physical_axis.raw_variable,
        "selection_axis": radial_axis.raw_variable,
        "status": status,
        "panels": panel_metadata,
        "transport_normalization": transport_normalization,
        "radial_normalization": radial_normalization,
        "weight_label": panel_weight.label,
        "weight_unit": panel_weight.unit,
        "cmap": cmap.name if cmap is not None else None,
        "norm": norm_metadata,
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _plot_conditioned_physical_pair_shells(
    interior: np.ndarray,
    axes: Sequence[DisplayAxis],
    radial_index: int,
    angular_index: int,
    x_index: int,
    y_index: int,
    sector: str,
    weight: WeightPresentation,
    weight_variable: str,
    *,
    product: str,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    _, _, apply_axes_style = _analysis_api(str(analysis_root))
    radial_axis = axes[radial_index]
    angular_axis = axes[angular_index]
    x_axis = axes[x_index]
    y_axis = axes[y_index]
    lower, upper = _angular_sector_bounds(config, sector)
    angular_weights = _angular_bin_selection_weights(angular_axis, lower, upper)
    angular_shape = [1] * interior.ndim
    angular_shape[angular_index] = angular_weights.size
    angular_selected = interior * angular_weights.reshape(angular_shape)

    scaled_panels = []
    panel_metadata = []
    panel_weight = weight
    transport_normalization = None
    radial_normalization = None
    radial_widths = np.diff(radial_axis.edges)
    radial_edges = radial_axis.edges
    shell_fraction = upper - lower
    for target_radius in FIXED_RADIAL_SHELL_TARGETS_KPC:
        selected_index = int(
            np.argmin(np.abs(np.log(radial_axis.centers / target_radius)))
        )
        radial_weights = np.zeros(radial_axis.centers.size, dtype=float)
        radial_weights[selected_index] = 1.0
        radial_shape = [1] * interior.ndim
        radial_shape[radial_index] = radial_weights.size
        marginal = _marginal_2d(
            angular_selected * radial_weights.reshape(radial_shape),
            x_index,
            y_index,
        )
        current_weight = weight
        if weight.style_key == "volume":
            shell_volume = (
                shell_fraction
                * (4.0 * math.pi / 3.0)
                * (
                    radial_edges[selected_index + 1] ** 3
                    - radial_edges[selected_index] ** 3
                )
            )
            marginal = marginal / shell_volume
            current_weight = _shell_volume_fraction_presentation()
            radial_normalization = "divided_by_angular_sector_shell_volume"
        elif _is_radial_transport(weight_variable):
            marginal = marginal / radial_widths[selected_index]
            current_weight = _crossing_presentation(weight, weight_variable)
            transport_normalization = "divided_by_radial_bin_width"
        if scaled_panels and panel_weight != current_weight:
            raise ValueError("Conditioned shell panels changed weight presentation")
        panel_weight = current_weight
        scaled = np.asarray(marginal, dtype=float) * panel_weight.factor
        scaled_panels.append(scaled)
        panel_metadata.append(
            {
                "target_radius_kpc": target_radius,
                "selected_radial_bin_index": selected_index,
                "selected_radial_center_kpc": float(
                    radial_axis.centers[selected_index]
                ),
                "selected_radial_edges_kpc": [
                    float(radial_edges[selected_index]),
                    float(radial_edges[selected_index + 1]),
                ],
                "minimum": float(np.nanmin(scaled)),
                "maximum": float(np.nanmax(scaled)),
                "sum": float(np.nansum(scaled)),
            }
        )

    stacked = np.stack(scaled_panels)
    status = "ok"
    norm_metadata = None
    norm: Optional[colors.Normalize] = None
    if np.any(stacked != 0.0):
        norm, plotted, norm_metadata = _weight_norm(stacked, panel_weight.signed)
        plotted_panels = [
            plotted[index] for index in range(len(FIXED_RADIAL_SHELL_TARGETS_KPC))
        ]
    else:
        status = "no_interior_support"
        plotted_panels = [np.ma.masked_all_like(values) for values in scaled_panels]

    fig, subplot_axes = plt.subplots(
        3, 2, figsize=(12.5, 13.0), squeeze=False, constrained_layout=True
    )
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    fig.suptitle(time_title, fontfamily="monospace")
    artist = None
    cmap = (
        _heatmap_cmap(config, panel_weight.style_key, norm)
        if norm is not None
        else None
    )
    for index, target_radius in enumerate(FIXED_RADIAL_SHELL_TARGETS_KPC):
        ax = subplot_axes.flat[index]
        ax.set_title(
            rf"$r \simeq {target_radius:g}\;\mathrm{{kpc}}$; "
            + ANGULAR_SECTOR_TITLES[sector],
            fontsize="small",
        )
        ax.set_xlabel(x_axis.label)
        ax.set_ylabel(y_axis.label)
        _apply_axis_scale(ax, x_axis, "x")
        _apply_axis_scale(ax, y_axis, "y")
        _decorate_axes(ax, source_context=source_context)
        apply_axes_style(ax, config["axes_style"])
        if norm is None:
            ax.set_facecolor(EMPTY_SUPPORT_COLOR)
            ax.text(
                0.5,
                0.5,
                "NO JOINT INTERIOR SUPPORT",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize="small",
                fontweight="bold",
            )
            continue
        artist = ax.pcolormesh(
            x_axis.edges,
            y_axis.edges,
            plotted_panels[index].T,
            shading="flat",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        _draw_reference_lines(ax, x_axis, "x", config)
        _draw_reference_lines(ax, y_axis, "y", config)
    if artist is not None:
        colorbar = fig.colorbar(
            artist, ax=list(subplot_axes.flat), fraction=0.025, pad=0.02
        )
        colorbar.set_label(_weight_axis_label(panel_weight))
        colorbar.ax.tick_params(
            which="both", direction=config["axes_style"]["tick_direction"]
        )
    paths = _save_figure(fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
    return {
        "kind": "conditioned_physical_pair_radial_shells",
        "paths": paths,
        "x_axis": x_axis.raw_variable,
        "y_axis": y_axis.raw_variable,
        "radial_axis": radial_axis.raw_variable,
        "angular_axis": angular_axis.raw_variable,
        "sector": sector,
        "abs_costheta_range": [lower, upper],
        "status": status,
        "panels": panel_metadata,
        "transport_normalization": transport_normalization,
        "radial_normalization": radial_normalization,
        "weight_label": panel_weight.label,
        "weight_unit": panel_weight.unit,
        "cmap": cmap.name if cmap is not None else None,
        "norm": norm_metadata,
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _thermokinematic_velocity_component(
    axis: DisplayAxis, velocity_sign: str
) -> Tuple[DisplayAxis, np.ndarray]:
    if axis.raw_variable != "velocity_r_km_s":
        raise ValueError("Thermokinematic velocity panel requires velocity_r_km_s")
    if velocity_sign == "positive":
        indices = np.flatnonzero(
            (axis.edges[:-1] >= 0.0) & (axis.edges[1:] > 0.0)
        )
        if indices.size == 0:
            raise ValueError("Thermokinematic velocity axis has no positive bins")
        edges = axis.edges[indices[0] : indices[-1] + 2]
        centers = axis.centers[indices]
        label = r"$v_r\;[\mathrm{km\,s^{-1}}]$"
    elif velocity_sign == "negative":
        indices = np.flatnonzero(
            (axis.edges[:-1] < 0.0) & (axis.edges[1:] <= 0.0)
        )[::-1]
        if indices.size == 0:
            raise ValueError("Thermokinematic velocity axis has no negative bins")
        edges = np.abs(axis.edges[indices[-1] : indices[0] + 2])[::-1]
        centers = np.abs(axis.centers[indices])
        label = r"$|v_r|\;[\mathrm{km\,s^{-1}}]$"
    else:
        raise ValueError(f"Unsupported thermokinematic velocity sign {velocity_sign!r}")
    return (
        replace(
            axis,
            raw_variable="absolute_velocity_r_km_s",
            label=label,
            scale="log",
            edges=np.asarray(edges, dtype=float),
            centers=np.asarray(centers, dtype=float),
            linthresh=1.0,
        ),
        indices,
    )


def _thermokinematic_velocity_signs(weight_variable: str) -> Tuple[str, ...]:
    """Return only velocity-sign panels that can carry this transport weight."""
    if weight_variable.endswith("_out"):
        return ("positive",)
    if weight_variable.endswith("_in") or weight_variable.endswith("_in_abs"):
        return ("negative",)
    return ("positive", "negative")


def _draw_thermokinematic_guides(
    ax: Any, velocity_axis: DisplayAxis, temperature_axis: DisplayAxis
) -> None:
    lower_velocity = max(
        THERMOKINEMATIC_VELOCITY_MIN_KM_S, float(velocity_axis.edges[0])
    )
    upper_velocity = min(
        THERMOKINEMATIC_VELOCITY_MAX_KM_S, float(velocity_axis.edges[-1])
    )
    velocities = np.geomspace(lower_velocity, upper_velocity, 256)
    temperature_lower = float(temperature_axis.edges[0])
    temperature_upper = float(temperature_axis.edges[-1])
    for mach in THERMOKINEMATIC_MACH_GUIDES:
        temperatures = (
            velocities / (mach * THERMOKINEMATIC_SOUND_SPEED_KM_S_PER_SQRT_K)
        ) ** 2
        valid = (temperatures >= temperature_lower) & (temperatures <= temperature_upper)
        if np.any(valid):
            ax.plot(
                velocities[valid],
                temperatures[valid],
                color="#252525",
                linestyle=":",
                linewidth=0.9,
                alpha=0.85,
            )
    for bernoulli_velocity in THERMOKINEMATIC_BERNOULLI_GUIDES_KM_S:
        temperatures = (
            (bernoulli_velocity**2 - velocities**2)
            * (THERMOKINEMATIC_GAMMA - 1.0)
            / (2.0 * THERMOKINEMATIC_SOUND_SPEED_KM_S_PER_SQRT_K**2)
        )
        valid = (temperatures >= temperature_lower) & (temperatures <= temperature_upper)
        if np.any(valid):
            ax.plot(
                velocities[valid],
                temperatures[valid],
                color="#b2182b",
                linestyle="--",
                linewidth=0.9,
                alpha=0.85,
            )
    ax.text(
        0.02,
        0.02,
        r"$\mathcal{M}=10^{-0.5},10^0,10^{0.5}$"
        "\n"
        r"$v_B=10^1,\ldots,10^{3.5}\;\mathrm{km\,s^{-1}}$",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize="xx-small",
        color="#252525",
        bbox={"facecolor": "white", "alpha": 0.76, "edgecolor": "#bdbdbd"},
    )


def _plot_thermokinematic_velocity_temperature_shells(
    interior: np.ndarray,
    axes: Sequence[DisplayAxis],
    radial_index: int,
    angular_index: int,
    temperature_index: int,
    velocity_index: int,
    sector: str,
    velocity_sign: str,
    weight: WeightPresentation,
    weight_variable: str,
    *,
    product: str,
    time_code: float,
    source_context: Mapping[str, Any],
    config: Mapping[str, Any],
    analysis_root: Path,
    output_base: Path,
    formats: Sequence[str],
    dpi: int,
    bbox_inches: str,
) -> Dict[str, Any]:
    _, _, apply_axes_style = _analysis_api(str(analysis_root))
    radial_axis = axes[radial_index]
    angular_axis = axes[angular_index]
    temperature_axis = axes[temperature_index]
    velocity_axis, velocity_indices = _thermokinematic_velocity_component(
        axes[velocity_index], velocity_sign
    )
    lower, upper = _angular_sector_bounds(config, sector)
    angular_weights = _angular_bin_selection_weights(angular_axis, lower, upper)
    angular_shape = [1] * interior.ndim
    angular_shape[angular_index] = angular_weights.size
    angular_selected = interior * angular_weights.reshape(angular_shape)

    scaled_panels = []
    panel_metadata = []
    panel_weight = weight
    transport_normalization = None
    radial_widths = np.diff(radial_axis.edges)
    radial_edges = radial_axis.edges
    for target_radius in FIXED_RADIAL_SHELL_TARGETS_KPC:
        selected_index = int(
            np.argmin(np.abs(np.log(radial_axis.centers / target_radius)))
        )
        radial_weights = np.zeros(radial_axis.centers.size, dtype=float)
        radial_weights[selected_index] = 1.0
        radial_shape = [1] * interior.ndim
        radial_shape[radial_index] = radial_weights.size
        marginal = _marginal_2d(
            angular_selected * radial_weights.reshape(radial_shape),
            velocity_index,
            temperature_index,
        )
        marginal = marginal[velocity_indices, :]
        current_weight = weight
        if _is_radial_transport(weight_variable):
            marginal = marginal / radial_widths[selected_index]
            current_weight = _crossing_presentation(weight, weight_variable)
            transport_normalization = "divided_by_radial_bin_width"
        if scaled_panels and panel_weight != current_weight:
            raise ValueError("Thermokinematic shell panels changed weight presentation")
        panel_weight = current_weight
        scaled = np.asarray(marginal, dtype=float) * panel_weight.factor
        scaled_panels.append(scaled)
        panel_metadata.append(
            {
                "target_radius_kpc": target_radius,
                "selected_radial_bin_index": selected_index,
                "selected_radial_center_kpc": float(
                    radial_axis.centers[selected_index]
                ),
                "selected_radial_edges_kpc": [
                    float(radial_edges[selected_index]),
                    float(radial_edges[selected_index + 1]),
                ],
                "minimum": float(np.nanmin(scaled)),
                "maximum": float(np.nanmax(scaled)),
                "sum": float(np.nansum(scaled)),
            }
        )

    stacked = np.stack(scaled_panels)
    status = "ok"
    norm_metadata = None
    norm: Optional[colors.Normalize] = None
    if np.any(stacked != 0.0):
        norm, plotted, norm_metadata = _weight_norm(stacked, panel_weight.signed)
        plotted_panels = [
            plotted[index] for index in range(len(FIXED_RADIAL_SHELL_TARGETS_KPC))
        ]
    else:
        status = "no_interior_support"
        plotted_panels = [np.ma.masked_all_like(values) for values in scaled_panels]

    fig, subplot_axes = plt.subplots(
        3, 2, figsize=(12.5, 13.0), squeeze=False, constrained_layout=True
    )
    time_title, time_since_feedback_myr = _time_title(time_code, config)
    sign_title = r"$v_r > 0$" if velocity_sign == "positive" else r"$v_r < 0$"
    fig.suptitle(f"{time_title}; {sign_title}", fontfamily="monospace")
    artist = None
    cmap = (
        _heatmap_cmap(config, panel_weight.style_key, norm)
        if norm is not None
        else None
    )
    for index, target_radius in enumerate(FIXED_RADIAL_SHELL_TARGETS_KPC):
        ax = subplot_axes.flat[index]
        ax.set_title(
            rf"$r \simeq {target_radius:g}\;\mathrm{{kpc}}$; "
            + ANGULAR_SECTOR_TITLES[sector],
            fontsize="small",
        )
        ax.set_xlabel(velocity_axis.label)
        ax.set_ylabel(temperature_axis.label)
        # The split velocity component retains the zero-adjacent source bin
        # edge for pcolormesh, but the displayed thermokinematic panel starts
        # at 1 km/s and therefore cannot use that zero edge as a log-axis
        # limit.
        ax.set_xscale("log")
        ax.set_xlim(
            THERMOKINEMATIC_VELOCITY_MIN_KM_S,
            THERMOKINEMATIC_VELOCITY_MAX_KM_S,
        )
        _apply_axis_scale(ax, temperature_axis, "y")
        _decorate_axes(ax, source_context=source_context)
        apply_axes_style(ax, config["axes_style"])
        if norm is None:
            ax.set_facecolor(EMPTY_SUPPORT_COLOR)
            ax.text(
                0.5,
                0.5,
                "NO JOINT INTERIOR SUPPORT",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize="small",
                fontweight="bold",
            )
            continue
        artist = ax.pcolormesh(
            velocity_axis.edges,
            temperature_axis.edges,
            plotted_panels[index].T,
            shading="flat",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        _draw_reference_lines(ax, temperature_axis, "y", config)
        _draw_thermokinematic_guides(ax, velocity_axis, temperature_axis)
    if artist is not None:
        colorbar = fig.colorbar(
            artist, ax=list(subplot_axes.flat), fraction=0.025, pad=0.02
        )
        colorbar.set_label(_weight_axis_label(panel_weight))
        colorbar.ax.tick_params(
            which="both", direction=config["axes_style"]["tick_direction"]
        )
    paths = _save_figure(fig, output_base, formats, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
    return {
        "kind": "thermokinematic_velocity_temperature_radial_shells",
        "paths": paths,
        "x_axis": velocity_axis.raw_variable,
        "y_axis": temperature_axis.raw_variable,
        "radial_axis": radial_axis.raw_variable,
        "angular_axis": angular_axis.raw_variable,
        "sector": sector,
        "velocity_sign": velocity_sign,
        "velocity_display_range_km_s": [
            THERMOKINEMATIC_VELOCITY_MIN_KM_S,
            THERMOKINEMATIC_VELOCITY_MAX_KM_S,
        ],
        "abs_costheta_range": [lower, upper],
        "status": status,
        "panels": panel_metadata,
        "transport_normalization": transport_normalization,
        "weight_label": panel_weight.label,
        "weight_unit": panel_weight.unit,
        "cmap": cmap.name if cmap is not None else None,
        "norm": norm_metadata,
        "guide_lines": {
            "mach": list(THERMOKINEMATIC_MACH_GUIDES),
            "bernoulli_km_s": list(THERMOKINEMATIC_BERNOULLI_GUIDES_KM_S),
            "sound_speed_km_s_per_sqrt_k": THERMOKINEMATIC_SOUND_SPEED_KM_S_PER_SQRT_K,
            "gamma": THERMOKINEMATIC_GAMMA,
        },
        "time_title": time_title,
        "time_since_feedback_myr": time_since_feedback_myr,
    }


def _source_fingerprint(
    product: str,
    payload_path: Path,
    header_path: Path,
    *,
    style_revision: str,
    style_path: Path,
    analysis_root: Path,
    formats: Sequence[str],
    dpi: int,
    profiles: bool,
    heatmaps: bool,
    source_context: Mapping[str, Any],
    snapshot_identity: Mapping[str, Any],
) -> Dict[str, Any]:
    payload_stat = payload_path.stat()
    header_stat = header_path.stat()
    return {
        "product": product,
        "payload_path": str(payload_path),
        "payload_size": payload_stat.st_size,
        "payload_mtime_ns": getattr(
            payload_stat, "st_mtime_ns", int(payload_stat.st_mtime * 1.0e9)
        ),
        "header_path": str(header_path),
        "header_size": header_stat.st_size,
        "header_mtime_ns": getattr(
            header_stat, "st_mtime_ns", int(header_stat.st_mtime * 1.0e9)
        ),
        "plotter_sha256": _sha256_file(Path(__file__)),
        "analysis_root": str(analysis_root),
        "analysis_reader_sha256": _sha256_file(
            analysis_root / "gotham_analysis/readers/pdf.py"
        ),
        "analysis_style_sha256": _sha256_file(
            analysis_root / "gotham_analysis/plotting/style.py"
        ),
        "style_path": str(style_path),
        "style_sha256": _sha256_file(style_path),
        "style_revision": style_revision,
        "matplotlib_version": matplotlib.__version__,
        "numpy_version": np.__version__,
        "pyyaml_version": yaml.__version__,
        "cmasher_version": getattr(cmasher, "__version__", None),
        "formats": list(formats),
        "dpi": dpi,
        "profiles": profiles,
        "heatmaps": heatmaps,
        "source_context": _json_ready(source_context),
        "snapshot_identity": _json_ready(snapshot_identity),
        "output_layout": "flat_scientific_movie_ready_v2",
    }


def _cached_product(
    metadata_path: Path, fingerprint: Mapping[str, Any]
) -> Optional[Dict[str, Any]]:
    try:
        with metadata_path.open("r", encoding="utf-8") as stream:
            metadata = json.load(stream)
    except (OSError, ValueError):
        return None
    if metadata.get("fingerprint") != _json_ready(fingerprint):
        return None
    paths = [
        Path(path)
        for plot in metadata.get("plots", [])
        for path in plot.get("paths", [])
    ]
    return metadata if paths and all(path.is_file() for path in paths) else None


def _remove_product_outputs(metadata_path: Path) -> None:
    try:
        with metadata_path.open("r", encoding="utf-8") as stream:
            metadata = json.load(stream)
    except (OSError, ValueError):
        metadata_path.unlink(missing_ok=True)
        return
    for plot in metadata.get("plots", []):
        for raw_path in plot.get("paths", []):
            Path(str(raw_path)).unlink(missing_ok=True)
    metadata_path.unlink(missing_ok=True)


def render_product(task: Mapping[str, Any]) -> Dict[str, Any]:
    product = str(task["product"])
    payload_path = Path(task["payload_path"])
    header_path = Path(task["header_path"])
    output_dir = Path(task["output_dir"])
    analysis_root = Path(task["analysis_root"])
    style_path = Path(task["style_path"])
    config = load_style_config(style_path)
    formats = tuple(task["formats"])
    profiles = bool(task["profiles"])
    heatmaps = bool(task["heatmaps"])
    source_context = dict(task["source_context"])
    snapshot_identity = dict(task["snapshot_identity"])
    image_dir = output_dir / "png"
    metadata_dir = output_dir / "metadata"
    metadata_path = metadata_dir / (
        _metadata_stem("product", snapshot_identity, product=product) + ".json"
    )
    fingerprint = _source_fingerprint(
        product,
        payload_path,
        header_path,
        style_revision=str(config["style_revision"]),
        style_path=style_path,
        analysis_root=analysis_root,
        formats=formats,
        dpi=int(task["dpi"]),
        profiles=profiles,
        heatmaps=heatmaps,
        source_context=source_context,
        snapshot_identity=snapshot_identity,
    )
    if task["skip_existing"]:
        cached = _cached_product(metadata_path, fingerprint)
        if cached is not None:
            cached["cached"] = True
            return cached
        _remove_product_outputs(metadata_path)
    elif task["overwrite"]:
        _remove_product_outputs(metadata_path)
    elif metadata_path.exists():
        raise FileExistsError(
            f"Plot metadata already exists: {metadata_path}; use --skip-existing "
            "or --overwrite"
        )

    read_pdf, _, _ = _analysis_api(str(analysis_root))
    snapshot = read_pdf(payload_path, header_path=header_path)
    if not math.isclose(
        float(snapshot["time"]),
        float(snapshot_identity["time_code"]),
        rel_tol=0.0,
        abs_tol=1.0e-10,
    ):
        raise ValueError(
            f"{product} payload time {snapshot['time']} does not match batch "
            f"time {snapshot_identity['time_code']}"
        )
    data = np.asarray(snapshot["pdf"], dtype=float)
    axes = [display_axis(info) for info in snapshot["header"]["dimensions"]]
    weight = weight_presentation(snapshot["header"])
    weight = _legacy_compatibility_weight(weight, product)
    if not np.all(np.isfinite(data)):
        raise ValueError(f"{product} contains non-finite histogram weights")
    if not weight.signed and np.any(data < 0.0):
        minimum = float(np.min(data))
        raise ValueError(
            f"{product} is a positive-weight product but contains negative "
            f"histogram weights (minimum {minimum:.17g})"
        )
    interior = data[tuple(slice(1, -1) for _ in range(data.ndim))]
    excluded = _excluded_summary(data, interior)
    dpi = int(task["dpi"])
    bbox_inches = str(config.get("output", {}).get("bbox_inches", "tight"))
    weight_variable = str(snapshot["header"].get("weight_variable", ""))
    is_transport = any(token in weight_variable for token in ("mdot", "edot", "ram"))

    plots: List[Dict[str, Any]] = []
    skipped_heatmaps: List[Dict[str, Any]] = []
    skipped_profiles: List[Dict[str, Any]] = []
    radial_index = _radial_axis_index(axes)
    angular_index = _angular_axis_index(axes)
    physical_indices = [
        index
        for index in range(len(axes))
        if index not in {radial_index, angular_index}
    ]
    conditioned_phase_recipe = (
        product.startswith("science_phase_")
        and not product.startswith("science_phase_opening_angle_")
        and radial_index is not None
        and angular_index is not None
        and len(axes) == 4
        and len(physical_indices) == 2
    )
    thermokinematic_recipe = (
        product.startswith("science_transport_thermokinematic_")
        and radial_index is not None
        and angular_index is not None
        and len(axes) == 4
        and len(physical_indices) == 2
    )
    custom_recipe = conditioned_phase_recipe or thermokinematic_recipe
    dedicated_2d_recipe = len(axes) == 2 and (
        product.startswith("science_transport_geometry_")
        or product.startswith("science_radial_mach_")
        or product.startswith("science_vertical_geometry_")
        or product.startswith("science_vertical_phase_")
        or product.startswith("science_angular_momentum_")
        or weight_variable
        in {
            "gas_metal_mdot_out",
            "gas_metal_mdot_in_abs",
            "total_metal_mdot_out",
            "total_metal_mdot_in_abs",
            "dust_mdot_out",
            "dust_mdot_in_abs",
        }
    )
    if profiles and custom_recipe:
        skipped_profiles.append(
            {
                "reason": (
                    "This four-dimensional product uses its dedicated "
                    "conditioned physical-pair plot recipe."
                )
            }
        )
    elif profiles and dedicated_2d_recipe:
        pass
    elif profiles:
        if radial_index is None or angular_index is None or not physical_indices:
            skipped_profiles.append(
                {
                    "reason": (
                        "Conditional mean/median profiles require radius, "
                        "polar angle, and at least one physical quantity axis."
                    )
                }
            )
        elif weight.signed:
            skipped_profiles.append(
                {
                    "reason": (
                        "Conditional means and medians are undefined for signed "
                        "histogram weights."
                    )
                }
            )
        else:
            for physical_index in physical_indices:
                name = _semantic_plot_name(
                    [
                        _axis_quantity_name(axes[physical_index].raw_variable),
                        "radius",
                        "costheta",
                    ],
                    axes,
                    {radial_index, angular_index, physical_index},
                    weight,
                )
                plots.append(
                    _plot_conditional_profile_summary(
                        interior,
                        axes,
                        radial_index,
                        angular_index,
                        physical_index,
                        time_code=float(snapshot["time"]),
                        source_context=source_context,
                        config=config,
                        analysis_root=analysis_root,
                        output_base=image_dir / _plot_stem(name, snapshot_identity),
                        formats=formats,
                        dpi=dpi,
                        bbox_inches=bbox_inches,
                    )
                )
    if heatmaps and custom_recipe:
        assert radial_index is not None
        assert angular_index is not None
        x_index, y_index = physical_indices
        retained_all = set(range(len(axes)))
        if thermokinematic_recipe:
            temperature_index = next(
                index
                for index in physical_indices
                if axes[index].raw_variable == "temperature_kelvin"
            )
            velocity_index = next(
                index
                for index in physical_indices
                if axes[index].raw_variable == "velocity_r_km_s"
            )
            for sector in ANGULAR_SECTOR_ORDER:
                for velocity_sign in _thermokinematic_velocity_signs(weight_variable):
                    name = _semantic_plot_name(
                        [
                            "radial_velocity",
                            "temperature",
                            "radius_shells",
                            sector,
                            f"{velocity_sign}_vr",
                        ],
                        axes,
                        retained_all,
                        weight,
                    )
                    plots.append(
                        _plot_thermokinematic_velocity_temperature_shells(
                            interior,
                            axes,
                            radial_index,
                            angular_index,
                            temperature_index,
                            velocity_index,
                            sector,
                            velocity_sign,
                            weight,
                            weight_variable,
                            product=product,
                            time_code=float(snapshot["time"]),
                            source_context=source_context,
                            config=config,
                            analysis_root=analysis_root,
                            output_base=image_dir / _plot_stem(name, snapshot_identity),
                            formats=formats,
                            dpi=dpi,
                            bbox_inches=bbox_inches,
                        )
                    )
        else:
            for sector in ANGULAR_SECTOR_ORDER:
                name = _semantic_plot_name(
                    [
                        _axis_quantity_name(axes[x_index].raw_variable),
                        _axis_quantity_name(axes[y_index].raw_variable),
                        "radius_shells",
                        sector,
                    ],
                    axes,
                    retained_all,
                    weight,
                )
                plots.append(
                    _plot_conditioned_physical_pair_shells(
                        interior,
                        axes,
                        radial_index,
                        angular_index,
                        x_index,
                        y_index,
                        sector,
                        weight,
                        weight_variable,
                        product=product,
                        time_code=float(snapshot["time"]),
                        source_context=source_context,
                        config=config,
                        analysis_root=analysis_root,
                        output_base=image_dir / _plot_stem(name, snapshot_identity),
                        formats=formats,
                        dpi=dpi,
                        bbox_inches=bbox_inches,
                    )
                )
        if thermokinematic_recipe:
            for physical_index, integrated_index in (
                (x_index, y_index),
                (y_index, x_index),
            ):
                name = _semantic_plot_name(
                    [
                        _axis_quantity_name(axes[physical_index].raw_variable),
                        "radius",
                        "theta_cuts",
                        _axis_quantity_name(
                            axes[integrated_index].raw_variable
                        )
                        + "_integrated",
                    ],
                    axes,
                    retained_all,
                    weight,
                )
                plots.append(
                    _plot_angular_sector_heatmaps(
                        interior,
                        axes,
                        radial_index,
                        physical_index,
                        angular_index,
                        weight,
                        weight_variable,
                        product=product,
                        time_code=float(snapshot["time"]),
                        source_context=source_context,
                        config=config,
                        analysis_root=analysis_root,
                        output_base=image_dir
                        / _plot_stem(name, snapshot_identity),
                        formats=formats,
                        dpi=dpi,
                        bbox_inches=bbox_inches,
                    )
                )
            name = _semantic_plot_name(
                [
                    "radius",
                    "costheta",
                    _axis_quantity_name(axes[x_index].raw_variable)
                    + "_integrated",
                    _axis_quantity_name(axes[y_index].raw_variable)
                    + "_integrated",
                ],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_heatmap(
                    _marginal_2d(interior, radial_index, angular_index),
                    axes[radial_index],
                    axes[angular_index],
                    weight,
                    weight_variable,
                    excluded,
                    product=product,
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )
    elif heatmaps:
        for first_index in range(len(axes)):
            for second_index in range(first_index + 1, len(axes)):
                x_index, y_index = _ordered_pair(first_index, second_index, axes)
                x_axis = axes[x_index]
                y_axis = axes[y_index]
                if _replaced_by_sectorized_radial_pair(
                    x_index,
                    y_index,
                    radial_index=radial_index,
                    angular_index=angular_index,
                ):
                    assert radial_index is not None
                    assert angular_index is not None
                    physical_index = (
                        y_index if x_index == angular_index else x_index
                    )
                    angular_axis = axes[angular_index]
                    physical_axis = axes[physical_index]
                    name = _semantic_plot_name(
                        [
                            _axis_quantity_name(physical_axis.raw_variable),
                            "costheta",
                            "radius_shells",
                        ],
                        axes,
                        {radial_index, angular_index, physical_index},
                        weight,
                    )
                    base = image_dir / _plot_stem(name, snapshot_identity)
                    plots.append(
                        _plot_fixed_radial_shell_heatmaps(
                            interior,
                            axes,
                            angular_index,
                            physical_index,
                            radial_index,
                            weight,
                            weight_variable,
                            product=product,
                            time_code=float(snapshot["time"]),
                            source_context=source_context,
                            config=config,
                            analysis_root=analysis_root,
                            output_base=base,
                            formats=formats,
                            dpi=dpi,
                            bbox_inches=bbox_inches,
                        )
                    )
                    continue
                if _redundant_state_geometry_pair(x_axis, y_axis, weight):
                    skipped_heatmaps.append(
                        {
                            "x_axis": x_axis.raw_variable,
                            "y_axis": y_axis.raw_variable,
                            "reason": (
                                "Mass- and volume-weighted radius--polar-angle "
                                "marginals reproduce weighted geometry rather "
                                "than the physical quantity carried by the "
                                "product."
                            ),
                        }
                    )
                    continue
                sectorized = _sectorized_radial_pair(
                    x_index,
                    y_index,
                    radial_index=radial_index,
                    angular_index=angular_index,
                )
                if sectorized:
                    physical_index = (
                        y_index if x_index == radial_index else x_index
                    )
                    name_components = [
                        _axis_quantity_name(axes[physical_index].raw_variable),
                        "radius",
                        "theta_cuts",
                    ]
                    retained_indices = {
                        radial_index,
                        angular_index,
                        physical_index,
                    }
                else:
                    name_components = [
                        _axis_quantity_name(x_axis.raw_variable),
                        _axis_quantity_name(y_axis.raw_variable),
                    ]
                    retained_indices = {x_index, y_index}
                name = _semantic_plot_name(
                    name_components,
                    axes,
                    retained_indices,
                    weight,
                )
                base = image_dir / _plot_stem(name, snapshot_identity)
                if sectorized:
                    assert angular_index is not None
                    plots.append(
                        _plot_angular_sector_heatmaps(
                            interior,
                            axes,
                            x_index,
                            y_index,
                            angular_index,
                            weight,
                            weight_variable,
                            product=product,
                            time_code=float(snapshot["time"]),
                            source_context=source_context,
                            config=config,
                            analysis_root=analysis_root,
                            output_base=base,
                            formats=formats,
                            dpi=dpi,
                            bbox_inches=bbox_inches,
                        )
                    )
                    continue
                plots.append(
                    _plot_heatmap(
                        _marginal_2d(interior, x_index, y_index),
                        x_axis,
                        y_axis,
                        weight,
                        weight_variable,
                        excluded,
                        product=product,
                        time_code=float(snapshot["time"]),
                        source_context=source_context,
                        config=config,
                        analysis_root=analysis_root,
                        output_base=base,
                        formats=formats,
                        dpi=dpi,
                        bbox_inches=bbox_inches,
                    )
                )

    if profiles and dedicated_2d_recipe:
        variable_indices = {
            axis.raw_variable: index for index, axis in enumerate(axes)
        }
        retained_all = set(range(len(axes)))

        if product.startswith("science_transport_geometry_"):
            radius_index = variable_indices["coord_r"]
            angle_index = variable_indices["coord_costheta"]
            angular_selections = []
            for sector in ANGULAR_SECTOR_ORDER:
                lower, upper = _angular_sector_bounds(config, sector)
                angular_selections.append(
                    (
                        sector,
                        ANGULAR_SECTOR_TITLES[sector],
                        _angular_bin_selection_weights(
                            axes[angle_index], lower, upper
                        ),
                    )
                )
            name = _semantic_plot_name(
                ["radius", "theta_sector_profiles"],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_selected_profiles_2d(
                    interior,
                    axes,
                    radius_index,
                    angle_index,
                    angular_selections,
                    weight,
                    weight_variable,
                    kind="angular_sector_radial_transport_profiles",
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir
                    / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )
            name = _semantic_plot_name(
                ["costheta", "fixed_radius_profiles"],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_fixed_primary_profiles_2d(
                    interior,
                    axes,
                    radius_index,
                    angle_index,
                    weight,
                    weight_variable,
                    kind="fixed_radius_angular_transport_profiles",
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir
                    / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )

        elif weight_variable in {
            "gas_metal_mdot_out",
            "gas_metal_mdot_in_abs",
            "total_metal_mdot_out",
            "total_metal_mdot_in_abs",
            "dust_mdot_out",
            "dust_mdot_in_abs",
        }:
            radius_index = variable_indices["coord_r"]
            temperature_index = variable_indices["temperature_kelvin"]
            temperature_selections = [
                (
                    "all",
                    "All temperatures",
                    np.ones(axes[temperature_index].centers.size),
                )
            ]
            for phase, limits in config.get("temperature_phases_k", {}).items():
                lower = None if limits[0] is None else float(limits[0])
                upper = None if limits[1] is None else float(limits[1])
                temperature_selections.append(
                    (
                        str(phase),
                        str(phase).replace("_", " ").title(),
                        _bin_range_selection_weights(
                            axes[temperature_index], lower, upper
                        ),
                    )
                )
            name = _semantic_plot_name(
                ["radius", "temperature_phase_profiles"],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_selected_profiles_2d(
                    interior,
                    axes,
                    radius_index,
                    temperature_index,
                    temperature_selections,
                    weight,
                    weight_variable,
                    kind="thermal_phase_radial_transport_profiles",
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir
                    / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )
            name = _semantic_plot_name(
                ["temperature", "fixed_radius_profiles"],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_fixed_primary_profiles_2d(
                    interior,
                    axes,
                    radius_index,
                    temperature_index,
                    weight,
                    weight_variable,
                    kind="fixed_radius_temperature_transport_profiles",
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir
                    / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )

        elif product.startswith("science_radial_mach_"):
            radius_index = variable_indices["coord_r"]
            mach_index = variable_indices["absolute_radial_mach"]
            mach_selections = [
                ("all", "All radial Mach numbers", np.ones(axes[mach_index].centers.size)),
                (
                    "subsonic",
                    r"$|{\cal M}_r| < 1$",
                    _bin_range_selection_weights(axes[mach_index], None, 1.0),
                ),
                (
                    "supersonic",
                    r"$|{\cal M}_r| \geq 1$",
                    _bin_range_selection_weights(axes[mach_index], 1.0, None),
                ),
            ]
            name = _semantic_plot_name(
                ["radius", "radial_mach_regime_profiles"],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_selected_profiles_2d(
                    interior,
                    axes,
                    radius_index,
                    mach_index,
                    mach_selections,
                    weight,
                    weight_variable,
                    kind="radial_mach_regime_profiles",
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir
                    / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )
            name = _semantic_plot_name(
                ["radial_mach", "fixed_radius_profiles"],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_fixed_primary_profiles_2d(
                    interior,
                    axes,
                    radius_index,
                    mach_index,
                    weight,
                    weight_variable,
                    kind="fixed_radius_radial_mach_profiles",
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir
                    / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )

        elif product.startswith("science_vertical_geometry_"):
            for variable, component in (
                ("cylindrical_radius", "cylindrical_radius_integrated_profile"),
                ("absolute_z", "absolute_z_integrated_profile"),
            ):
                profile_index = variable_indices[variable]
                name = _semantic_plot_name(
                    [component], axes, retained_all, weight
                )
                plots.append(
                    _plot_integrated_profile_1d(
                        interior,
                        axes,
                        profile_index,
                        weight,
                        kind="vertical_geometry_integrated_profile",
                        time_code=float(snapshot["time"]),
                        source_context=source_context,
                        config=config,
                        analysis_root=analysis_root,
                        output_base=image_dir
                        / _plot_stem(name, snapshot_identity),
                        formats=formats,
                        dpi=dpi,
                        bbox_inches=bbox_inches,
                    )
                )

        elif product.startswith("science_vertical_phase_"):
            height_index = variable_indices["absolute_z"]
            temperature_index = variable_indices["temperature_kelvin"]
            temperature_selections = []
            for phase, limits in config.get("temperature_phases_k", {}).items():
                lower = None if limits[0] is None else float(limits[0])
                upper = None if limits[1] is None else float(limits[1])
                temperature_selections.append(
                    (
                        str(phase),
                        str(phase).replace("_", " ").title(),
                        _bin_range_selection_weights(
                            axes[temperature_index], lower, upper
                        ),
                    )
                )
            name = _semantic_plot_name(
                ["absolute_z", "temperature_phase_profiles"],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_selected_profiles_2d(
                    interior,
                    axes,
                    height_index,
                    temperature_index,
                    temperature_selections,
                    weight,
                    weight_variable,
                    kind="vertical_thermal_phase_profiles",
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir
                    / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )

        elif product.startswith("science_angular_momentum_"):
            radius_index = variable_indices["coord_r"]
            angular_momentum_index = variable_indices[
                "specific_angular_momentum_z_kpc_km_s"
            ]
            name = _semantic_plot_name(
                ["radius", "angular_momentum_integrated_profile"],
                axes,
                retained_all,
                weight,
            )
            plots.append(
                _plot_selected_profiles_2d(
                    interior,
                    axes,
                    radius_index,
                    angular_momentum_index,
                    [
                        (
                            "all",
                            "All angular momentum",
                            np.ones(
                                axes[angular_momentum_index].centers.size
                            ),
                        )
                    ],
                    weight,
                    weight_variable,
                    kind="radial_angular_momentum_transport_profile",
                    time_code=float(snapshot["time"]),
                    source_context=source_context,
                    config=config,
                    analysis_root=analysis_root,
                    output_base=image_dir
                    / _plot_stem(name, snapshot_identity),
                    formats=formats,
                    dpi=dpi,
                    bbox_inches=bbox_inches,
                )
            )

    metadata = {
        "schema_version": 1,
        "kind": "rebuilt_pdf_plot_product",
        "product": product,
        "cached": False,
        "fingerprint": fingerprint,
        "source_payload": str(payload_path),
        "source_header": str(header_path),
        "time_code": float(snapshot["time"]),
        "snapshot_identity": snapshot_identity,
        "metadata_path": str(metadata_path),
        "source_context": source_context,
        "variables": [axis.raw_variable for axis in axes],
        "shape_with_overflow": list(data.shape),
        "interior_shape": list(interior.shape),
        "distribution_kind": "bin-integrated weighted histogram",
        "marginalization": (
            "All axes are restricted to interior bins before summation. "
            "Four-panel radial marginals use fractional polar-bin overlap for "
            "all-angle, polar, midplane, and intermediate sectors. "
            "Six-panel angular marginals select the nearest radial PDF bin to "
            "2, 4, 8, 16, 32, and 64 kpc. "
            "Volume-weighted radial panels are divided by the corresponding "
            "angular sector's analytic spherical-shell volume; fixed-shell "
            "angular panels are divided by each radial/angular wedge volume. "
            "Fixed-radius profile panels are divided by the displayed "
            "distribution bin width, using log10-bin width for log-scaled "
            "distribution axes. "
            "Other panels are not probability-density-normalized distributions."
        ),
        "axes": [_axis_metadata(axis) for axis in axes],
        "weight": {
            "header_weight": snapshot["header"].get("weight"),
            "header_weight_variable": weight_variable or None,
            "style_key": weight.style_key,
            "label": weight.label,
            "unit": weight.unit,
            "factor": weight.factor,
            "signed": weight.signed,
            "is_crossing_rate": False if is_transport else None,
            "transport_note": (
                "Transport weights are raw volume-integrated moments with an "
                "extra length dimension. Plots retaining radius divide by "
                "radial-bin width and are crossing-rate diagnostics; other "
                "transport plots remain raw integrated moments."
                if is_transport
                else None
            ),
            "radial_volume_note": (
                "Volume-weighted panels retaining coord_r divide each radial "
                "bin by its selected angular sector's fraction of "
                "4*pi/3*(r_outer^3-r_inner^3), so colors show the filling "
                "fraction within that sector. The all-angle radial profile "
                "uses the full shell volume. Fixed-shell angular panels divide "
                "by each radial/angular wedge volume."
                if weight.style_key == "volume"
                and any(axis.raw_variable == "coord_r" for axis in axes)
                else None
            ),
        },
        "excluded_bins": excluded,
        "overflow_by_axis": _overflow_by_axis(data, axes),
        "skipped_heatmaps": skipped_heatmaps,
        "skipped_profiles": skipped_profiles,
        "plots": plots,
        "plot_count": len(plots),
    }
    _write_json_atomic(metadata_path, metadata)
    return metadata


def load_rebuild_manifest(input_dir: Path) -> Dict[str, Any]:
    manifest_path = input_dir / "rebuild_manifest.json"
    if not manifest_path.is_file():
        raise ValueError(
            f"Completed reducer manifest is missing: {manifest_path}. "
            "Refusing to plot a potentially partial reconstruction."
        )
    with manifest_path.open("r", encoding="utf-8") as stream:
        manifest = json.load(stream)
    if not isinstance(manifest, dict):
        raise ValueError(
            f"Reducer manifest must contain a JSON object: {manifest_path}"
        )
    return manifest


def discover_products(input_dir: Path, patterns: Sequence[str]) -> List[Dict[str, str]]:
    manifest = load_rebuild_manifest(input_dir)
    product_ids = [
        str(product["id"])
        for product in manifest.get("products", [])
        if isinstance(product, dict) and product.get("id")
    ]
    product_ids = sorted(dict.fromkeys(product_ids))
    if patterns:
        product_ids = [
            product
            for product in product_ids
            if any(fnmatch.fnmatch(product, pattern) for pattern in patterns)
        ]
    discovered = []
    for product in product_ids:
        product_dir = input_dir / product
        header = product_dir / "gotham.header.pdf"
        payloads = (
            sorted(
                path
                for path in product_dir.iterdir()
                if path.is_file() and PAYLOAD_RE.match(path.name)
            )
            if product_dir.is_dir()
            else []
        )
        if not header.is_file() or len(payloads) != 1:
            raise ValueError(
                f"{product_dir} must contain one gotham.#####.pdf payload and gotham.header.pdf"
            )
        discovered.append(
            {
                "product": product,
                "payload_path": str(payloads[0]),
                "header_path": str(header),
            }
        )
    if not discovered:
        raise ValueError(f"No reconstructed PDF products selected beneath {input_dir}")
    return discovered


def _payload_time_code(payload_path: Path) -> float:
    with payload_path.open("rb") as stream:
        values = np.fromfile(stream, dtype=np.float64, count=1)
    if values.size != 1 or not np.isfinite(values[0]):
        raise ValueError(f"Cannot read finite snapshot time from {payload_path}")
    return float(values[0])


def _gallery_html(metadata: Mapping[str, Any], output_dir: Path) -> str:
    sections = []
    for product in metadata["products"]:
        product_name = html.escape(str(product["product"]))
        images = []
        for plot in product.get("plots", []):
            pngs = [
                path for path in plot.get("paths", []) if str(path).endswith(".png")
            ]
            if not pngs:
                continue
            relative = Path(pngs[0]).relative_to(output_dir)
            caption = plot["kind"]
            if plot["kind"] == "two_dimensional_marginal":
                caption = f"{plot['x_axis']} vs {plot['y_axis']}"
                if plot.get("radial_normalization") == (
                    "divided_by_angular_sector_shell_volume"
                ):
                    caption += " | shell-volume fraction"
                elif plot.get("transport_normalization"):
                    caption += " | radial-bin-width normalized"
                if plot.get("status") != "ok":
                    caption += f" | {plot.get('status')}"
            elif plot["kind"] == "angular_sector_radial_marginals":
                caption = (
                    f"{plot['x_axis']} vs {plot['y_axis']} | four angular sectors"
                )
            elif plot["kind"] == "fixed_radial_shell_angular_marginals":
                caption = (
                    f"{plot['x_axis']} vs {plot['y_axis']} | six radial shells"
                )
            elif plot["kind"] == "conditional_mean_median_profiles":
                caption = (
                    f"mean/median {plot['physical_axis']} vs radius and polar angle"
                )
            images.append(
                '<figure><a href="{0}"><img src="{0}" loading="lazy"></a>'
                "<figcaption>{1}</figcaption></figure>".format(
                    html.escape(str(relative)), html.escape(caption)
                )
            )
        sections.append(
            '<section><h2>{}</h2><div class="grid">{}</div></section>'.format(
                product_name, "".join(images)
            )
        )
    source_context = metadata.get("source_context", {})
    source_notice = ""
    if source_context.get("annotation"):
        source_notice = '<p class="warning">{}</p>'.format(
            html.escape(str(source_context["annotation"]))
        )
    return """<!doctype html>
<html><head><meta charset="utf-8"><title>GOTHAM rebuilt PDF plots</title>
<style>
body {{ font-family: sans-serif; margin: 2rem; background: #f5f5f5; color: #222; }}
.warning {{ color: #a50f15; font-weight: bold; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit,minmax(320px,1fr)); gap: 1rem; }}
figure {{ background: white; margin: 0; padding: .6rem; border: 1px solid #ccc; }}
img {{ width: 100%; height: auto; }} figcaption {{ margin-top: .4rem; font-family: monospace; }}
</style></head><body><h1>GOTHAM rebuilt PDF plots</h1>{}{}</body></html>
""".format(
        source_notice, "".join(sections)
    )


def _write_gallery(output_dir: Path, metadata: Mapping[str, Any]) -> None:
    path = output_dir / "index.html"
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(_gallery_html(metadata, output_dir), encoding="utf-8")
    temporary.replace(path)


def _default_output_dir(input_dir: Path) -> Path:
    return input_dir.with_name(input_dir.name + "_plots")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_dir", type=Path, help="One completed reducer output directory"
    )
    parser.add_argument(
        "--output-dir", type=Path, help="Separate plot tree; defaults to INPUT_plots"
    )
    parser.add_argument("--analysis-root", type=Path, default=DEFAULT_ANALYSIS_ROOT)
    parser.add_argument("--style-config", type=Path, default=DEFAULT_STYLE_CONFIG)
    parser.add_argument(
        "--render-profile",
        choices=("preview", "production", "paper"),
        default="preview",
    )
    parser.add_argument(
        "--format", dest="formats", action="append", choices=("png", "pdf")
    )
    parser.add_argument(
        "--dpi", type=int, help="Override the configured render-profile DPI"
    )
    parser.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument(
        "--product",
        dest="patterns",
        action="append",
        help="Product id or shell-style pattern; repeatable. Default: all products.",
    )
    existing = parser.add_mutually_exclusive_group()
    existing.add_argument("--skip-existing", action="store_true")
    existing.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-heatmaps", action="store_true")
    parser.add_argument("--no-profiles", action="store_true")
    parser.add_argument("--no-gallery", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    input_dir = args.input_dir.resolve()
    output_dir = (args.output_dir or _default_output_dir(input_dir)).resolve()
    analysis_root = args.analysis_root.resolve()
    style_path = args.style_config.resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {input_dir}")
    if (
        output_dir == input_dir
        or input_dir in output_dir.parents
        or output_dir in input_dir.parents
    ):
        raise SystemExit(
            "Refusing nested input/output trees; plots must use a separate sibling tree"
        )
    if args.no_heatmaps and args.no_profiles:
        raise SystemExit("At least one of heatmaps or profiles must be enabled")
    if args.workers <= 0:
        raise SystemExit("--workers must be positive")

    config = load_style_config(style_path)
    profile = config["render_profiles"][args.render_profile]
    formats = args.formats or list(profile["formats"])
    dpi = int(args.dpi or profile["dpi"])
    manifest = load_rebuild_manifest(input_dir)
    source_context = _source_context(manifest)
    products = discover_products(input_dir, args.patterns or [])
    time_code = _payload_time_code(Path(products[0]["payload_path"]))
    snapshot_identity = _snapshot_identity(
        input_dir, manifest, time_code=time_code, config=config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    _lock_stream = _acquire_output_lock(output_dir)

    common = {
        "output_dir": str(output_dir),
        "analysis_root": str(analysis_root),
        "style_path": str(style_path),
        "formats": formats,
        "dpi": dpi,
        "profiles": not args.no_profiles,
        "heatmaps": not args.no_heatmaps,
        "skip_existing": args.skip_existing,
        "overwrite": args.overwrite,
        "source_context": source_context,
        "snapshot_identity": snapshot_identity,
    }
    tasks = [dict(common, **product) for product in products]
    results: List[Dict[str, Any]] = []
    failures: List[Dict[str, str]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_tasks = {executor.submit(render_product, task): task for task in tasks}
        for future in concurrent.futures.as_completed(future_tasks):
            task = future_tasks[future]
            try:
                result = future.result()
            except Exception as exc:
                failures.append(
                    {
                        "product": str(task["product"]),
                        "error": str(exc),
                        "traceback": traceback.format_exc(),
                    }
                )
                print(f"FAIL {task['product']}: {exc}", file=sys.stderr)
                continue
            results.append(result)
            marker = "CACHED" if result.get("cached") else "DONE"
            print(f"{marker} {result['product']}: {result['plot_count']} plots")

    results.sort(key=lambda value: value["product"])
    metadata = {
        "schema_version": 1,
        "kind": "rebuilt_pdf_plot_batch",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "analysis_root": str(analysis_root),
        "style_config": str(style_path),
        "style_revision": config["style_revision"],
        "render_profile": args.render_profile,
        "formats": formats,
        "dpi": dpi,
        "workers": args.workers,
        "source_context": source_context,
        "snapshot_identity": snapshot_identity,
        "requested_product_patterns": args.patterns or [],
        "product_count": len(results),
        "failed_product_count": len(failures),
        "plot_count": sum(int(result["plot_count"]) for result in results),
        "distribution_note": (
            "Panels are bin-integrated weighted histograms with a recorded "
            "absolute-weight visual-support floor. All axes are restricted "
            "to interior bins before marginalization. Volume-weighted panels "
            "retaining radius are divided by the selected angular sector's "
            "analytic spherical-shell volume. "
            "Radial transport panels retaining radius are divided by radial-bin "
            "width. Fixed-radius profile panels are divided by their displayed "
            "distribution-bin width, using log10-bin width for log-scaled "
            "distribution axes; other transport panels remain raw integrated "
            "moments."
        ),
        "products": results,
        "failures": failures,
    }
    batch_metadata_path = output_dir / "metadata" / (
        _metadata_stem("batch", snapshot_identity) + ".json"
    )
    metadata["metadata_path"] = str(batch_metadata_path)
    _write_json_atomic(batch_metadata_path, metadata)
    _write_json_atomic(output_dir / "metadata.json", metadata)
    if not args.no_gallery:
        _write_gallery(output_dir, metadata)
    else:
        (output_dir / "index.html").unlink(missing_ok=True)
    print(
        f"Wrote {metadata['plot_count']} plot groups for {metadata['product_count']} "
        f"products to {output_dir}"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
