"""Independent pure-Python reference model for the standalone GOTHAM reducer."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

from synthetic_gotham import BlockSpec, GAMMA_DEFAULT, cell_center, iter_cells


@dataclass(frozen=True)
class Axis:
    variable: str
    nbin: int
    minimum: float
    maximum: float
    scale: str = "linear"
    linthresh: float = 1.0


@dataclass(frozen=True)
class Product:
    identifier: str
    axes: Tuple[Axis, ...]
    weight: str

    @property
    def shape(self) -> Tuple[int, ...]:
        return tuple(axis.nbin + 2 for axis in self.axes)

    @property
    def strides(self) -> Tuple[int, ...]:
        result = []
        stride = 1
        for axis in reversed(self.axes):
            result.append(stride)
            stride *= axis.nbin + 2
        return tuple(reversed(result))

    @property
    def total_bins(self) -> int:
        result = 1
        for extent in self.shape:
            result *= extent
        return result


def linear(variable: str, nbin: int, minimum: float, maximum: float) -> Axis:
    return Axis(variable, nbin, minimum, maximum)


def log_axis(variable: str, nbin: int, minimum: float, maximum: float) -> Axis:
    return Axis(variable, nbin, minimum, maximum, "log")


def symlog_axis(
    variable: str, nbin: int, minimum: float, maximum: float, linthresh: float
) -> Axis:
    return Axis(variable, nbin, minimum, maximum, "symlog", linthresh)


R32_WIDE = log_axis("coord_r", 32, 1.0e-1, 190.0)
R128_WIDE = log_axis("coord_r", 128, 1.0e-1, 190.0)
R32_CORE = log_axis("coord_r", 32, 1.0e-1, 20.0)
R4_CORE = log_axis("coord_r", 4, 1.0e-1, 20.0)
ABSCT16 = linear("coord_abscostheta", 16, 0.0, 1.0)
ABSCT4 = linear("coord_abscostheta", 4, 0.0, 1.0)
COST16 = linear("coord_costheta", 16, 0.0, 1.0)
TEMP128 = log_axis("temperature", 128, 1.0e-3, 7.0e4)
RHO128 = log_axis("hydro_w_d", 128, 1.0e-5, 1.0e5)
EINT128 = log_axis("hydro_w_e", 128, 5.0e-5, 5.0e3)
VR256 = linear("vel_sph_r", 256, -100.0, 500.0)
VT128 = linear("vel_sph_theta", 128, -25.0, 25.0)
VP128 = linear("vel_sph_phi", 128, -25.0, 25.0)
VR128_POSITIVE = log_axis("vel_sph_r", 128, 1.0e-1, 1.0e3)
VR64_POSITIVE = log_axis("vel_sph_r", 64, 1.0e-1, 1.0e3)
S0_128 = log_axis("hydro_w_s_00", 128, 1.0e-5, 0.3)
S1_128 = log_axis("hydro_w_s_01", 128, 1.0e-5, 1.0)
S2_128 = log_axis("hydro_w_s_02", 128, 1.0e-5, 1.0)

ORIGINAL_PRODUCTS = (
    Product("output23", (R32_WIDE, ABSCT16, TEMP128, RHO128), "mass"),
    Product("output24", (R128_WIDE, ABSCT16, EINT128), "volume"),
    Product("output25", (R128_WIDE, ABSCT16, VR256), "mass"),
    Product("output26", (R128_WIDE, ABSCT16, VP128), "mass"),
    Product("output27", (R128_WIDE, ABSCT16, VT128), "mass"),
    Product("output28", (R128_WIDE, ABSCT16), "mdot_sph"),
    Product("output29", (ABSCT16, R128_WIDE), "edot_sph"),
    Product("output30", (R32_CORE, ABSCT16, TEMP128, VR128_POSITIVE), "mdot_sph"),
    Product("output31", (COST16, R32_CORE, TEMP128, VR128_POSITIVE), "edot_sph"),
    Product("output32", (R32_CORE, ABSCT16, S0_128, VR64_POSITIVE), "mass"),
    Product("output33", (R32_CORE, ABSCT16, S1_128), "mass"),
    Product("output34", (R32_CORE, ABSCT16, S2_128), "mass"),
    Product("output35", (R4_CORE, ABSCT4, S1_128, S0_128), "mass"),
    Product("output36", (R4_CORE, ABSCT4, S2_128, S0_128), "mass"),
)

R_FULL = log_axis("coord_r", 361, 6.4e-2, 262.144)
ABS_MU = linear("coord_abscostheta", 64, 0.0, 1.0)
SIGNED_MU = linear("coord_costheta", 128, -1.0, 1.0)
R_PHASE = log_axis("coord_r", 64, 6.4e-2, 128.0)
ABS_MU_PHASE = linear("coord_abscostheta", 16, 0.0, 1.0)
Z_FULL = log_axis("absolute_z", 192, 2.0e-3, 210.0)
RCYL_FULL = log_axis("cylindrical_radius", 192, 2.0e-3, 290.0)
PRESSURE_KB = log_axis("pressure_over_kb_k_cm3", 256, 1.0e-6, 1.0e10)
NH_PHASE = log_axis("hydrogen_number_density_cm3", 96, 1.0e-10, 1.0e6)
TEMP_PHASE = log_axis("temperature_kelvin", 96, 1.0, 1.0e9)
TEMP_TRANSPORT = log_axis("temperature_kelvin", 80, 1.0, 1.0e9)
VR_TRANSPORT = symlog_axis("velocity_r_km_s", 80, -1.0e4, 1.0e4, 5.0)
VTHETA = symlog_axis("velocity_theta_km_s", 192, -2.5e3, 2.5e3, 5.0)
TEMP_K = log_axis("temperature_kelvin", 192, 1.0, 1.0e9)
COOLING_RATE = symlog_axis(
    "cooling_rate_erg_s_cm3", 192, -1.0e-14, 1.0e-14, 1.0e-30
)
COOLING_TIME = symlog_axis("cooling_time_myr", 192, -1.0e8, 1.0e8, 1.0e-1)
JZ = symlog_axis(
    "specific_angular_momentum_z_kpc_km_s", 192, -1.0e6, 1.0e6, 10.0
)

# Representative products used for executable-level numerical tests. The full
# 69-product inventory is validated separately through --list-products so unit
# tests do not emit multi-gigabyte dense payloads.
SCIENCE_PRODUCTS = (
    Product(
        "science_r_theta_pressure_volume",
        (R_FULL, ABS_MU, PRESSURE_KB),
        "volume",
    ),
    Product(
        "science_r_theta_vtheta_mass",
        (R_FULL, ABS_MU, VTHETA),
        "mass",
    ),
    Product(
        "science_r_theta_edot_cool_volume",
        (R_FULL, ABS_MU, COOLING_RATE),
        "volume",
    ),
    Product(
        "science_r_theta_temperature_edot_cool",
        (R_FULL, ABS_MU, TEMP_K),
        "edot_cool",
    ),
    Product(
        "science_r_theta_tcool_volume",
        (R_FULL, ABS_MU, COOLING_TIME),
        "volume",
    ),
    Product(
        "science_phase_density_temperature_volume",
        (R_PHASE, ABS_MU_PHASE, NH_PHASE, TEMP_PHASE),
        "volume",
    ),
    Product(
        "science_phase_density_temperature_edot_cool",
        (R_PHASE, ABS_MU_PHASE, NH_PHASE, TEMP_PHASE),
        "edot_cool",
    ),
    Product(
        "science_transport_thermokinematic_edot_in_abs",
        (R_PHASE, ABS_MU_PHASE, TEMP_TRANSPORT, VR_TRANSPORT),
        "edot_sph_in_abs",
    ),
    Product(
        "science_transport_thermokinematic_edot_cool",
        (R_PHASE, ABS_MU_PHASE, TEMP_TRANSPORT, VR_TRANSPORT),
        "edot_cool",
    ),
    Product(
        "science_gas_metal_mdot_in_abs",
        (R_FULL, TEMP_K),
        "gas_metal_mdot_in_abs",
    ),
    Product(
        "science_vertical_geometry_ram_out",
        (RCYL_FULL, Z_FULL),
        "vertical_ram_out",
    ),
    Product(
        "science_angular_momentum_mdot_out",
        (R_FULL, JZ),
        "mdot_sph_out",
    ),
)

SCIENCE_PRODUCT_COUNT = 69
SCIENCE_TOTAL_BINS = 305_811_772

LENGTH_CGS = 3.0856775809623245e21
MASS_CGS = 3.036951775493658e39
TIME_CGS = 3.15576e15
MYR_CGS = 3.15576e13
MU = 0.62
K_BOLTZMANN_CGS = 1.380649e-16
PROTON_MASS_CGS = 1.67262192369e-24
VELOCITY_KM_S = LENGTH_CGS / TIME_CGS / 1.0e5
DENSITY_CGS = MASS_CGS / LENGTH_CGS**3
PRESSURE_CGS = MASS_CGS / (LENGTH_CGS * TIME_CGS**2)
TEMPERATURE_K_PER_RATIO = (
    (GAMMA_DEFAULT - 1.0)
    * PRESSURE_CGS
    / DENSITY_CGS
    * MU
    * PROTON_MASS_CGS
    / K_BOLTZMANN_CGS
)
PRESSURE_OVER_KB_PER_EINT = (
    (GAMMA_DEFAULT - 1.0) * PRESSURE_CGS / K_BOLTZMANN_CGS
)
NH_PER_DENSITY = 0.0463500
COOLING_HRATE = 2.0e-26
COOLING_HSCALE_NORM = 1347.6541799847946
COOLING_HSCALE_RADIUS = 0.9402361259893334
COOLING_HSCALE_HEIGHT = 0.04950000000848683
COOLING_SOLAR_METALLICITY = 0.02
FLOAT32_TINY = 1.1754943508222875e-38


def symlog_forward(value: float, linthresh: float) -> float:
    absolute = abs(value)
    transformed = absolute / linthresh if absolute <= linthresh else 1.0 + math.log10(absolute / linthresh)
    return math.copysign(transformed, value)


def transform(value: float, axis: Axis) -> float:
    if axis.scale == "log":
        return math.log10(value)
    if axis.scale == "symlog":
        return symlog_forward(value, axis.linthresh)
    return value


def bin_index(value: float, axis: Axis) -> int:
    if value < axis.minimum:
        return 0
    if not value < axis.maximum:
        return axis.nbin + 1
    position = (transform(value, axis) - transform(axis.minimum, axis)) * axis.nbin / (
        transform(axis.maximum, axis) - transform(axis.minimum, axis)
    )
    if position < 0.0:
        return 0
    if position >= axis.nbin:
        return axis.nbin
    return int(position) + 1


@lru_cache(maxsize=1)
def cooling_tables():
    """Read the exact table payload that the standalone reducer includes."""

    path = Path(__file__).resolve().parents[3] / "src" / "srcterms" / "cooling_tables.hpp"
    text = path.read_text(encoding="utf-8")

    def array(name: str) -> Tuple[float, ...]:
        match = re.search(
            rf"constexpr Real {re.escape(name)}\[.*?\] = \{{(.*?)\}};",
            text,
            flags=re.DOTALL,
        )
        if match is None:
            raise RuntimeError(f"Could not find {name} in {path}")
        return tuple(
            float(value)
            for value in re.findall(
                r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", match.group(1)
            )
        )

    return {
        "temperature": array("Tbins_ARR"),
        "hydrogen": array("nHbins_ARR"),
        "metal_pie": array("Metal_Cooling_ARR"),
        "hhe_pie": array("H_He_Cooling_ARR"),
        "metal_cie": array("Metal_Cooling_CIE_ARR"),
        "hhe_cie": array("H_He_Cooling_CIE_ARR"),
    }


def cooling_lambdas(temperature: float, n_h: float, metallicity: float) -> Tuple[float, float]:
    tables = cooling_tables()
    temperature_bins = tables["temperature"]
    hydrogen_bins = tables["hydrogen"]
    log_temperature = math.log10(temperature)
    log_n_h = math.log10(n_h)
    temperature_floor = 10.0 ** temperature_bins[0]
    temperature_ceiling = 10.0 ** temperature_bins[-1]
    hydrogen_floor = 10.0 ** hydrogen_bins[0]
    hydrogen_ceiling = 10.0 ** hydrogen_bins[-1]
    low_temperature = 1.0 if temperature < temperature_floor else 0.0
    temperature_in_table = 1.0 if temperature_floor <= temperature <= temperature_ceiling else 0.0
    hydrogen_in_table = 1.0 if hydrogen_floor <= n_h <= hydrogen_ceiling else 0.0
    pie_in_table = temperature_in_table * hydrogen_in_table

    temperature_index = 0
    hydrogen_index = 0
    while temperature_index < len(temperature_bins) - 2 and temperature_bins[
        temperature_index + 1
    ] < log_temperature:
        temperature_index += 1
    while hydrogen_index < len(hydrogen_bins) - 2 and hydrogen_bins[
        hydrogen_index + 1
    ] < log_n_h:
        hydrogen_index += 1

    temperature_weight = (log_temperature - temperature_bins[temperature_index]) / (
        temperature_bins[temperature_index + 1] - temperature_bins[temperature_index]
    )
    inverse_temperature_weight = 1.0 - temperature_weight
    hydrogen_weight = (log_n_h - hydrogen_bins[hydrogen_index]) / (
        hydrogen_bins[hydrogen_index + 1] - hydrogen_bins[hydrogen_index]
    )
    inverse_hydrogen_weight = 1.0 - hydrogen_weight
    n_h_count = len(hydrogen_bins)
    pie00 = temperature_index * n_h_count + hydrogen_index
    pie10 = (temperature_index + 1) * n_h_count + hydrogen_index
    pie01 = temperature_index * n_h_count + hydrogen_index + 1
    pie11 = (temperature_index + 1) * n_h_count + hydrogen_index + 1

    primordial_pie = (
        inverse_temperature_weight * inverse_hydrogen_weight * tables["hhe_pie"][pie00]
        + temperature_weight * inverse_hydrogen_weight * tables["hhe_pie"][pie10]
        + inverse_temperature_weight * hydrogen_weight * tables["hhe_pie"][pie01]
        + temperature_weight * hydrogen_weight * tables["hhe_pie"][pie11]
    )
    metal_pie = (
        inverse_temperature_weight * inverse_hydrogen_weight * tables["metal_pie"][pie00]
        + temperature_weight * inverse_hydrogen_weight * tables["metal_pie"][pie10]
        + inverse_temperature_weight * hydrogen_weight * tables["metal_pie"][pie01]
        + temperature_weight * hydrogen_weight * tables["metal_pie"][pie11]
    )
    lambda_pie = pie_in_table * (metallicity * metal_pie + primordial_pie)

    primordial_cie = tables["hhe_cie"][temperature_index] + temperature_weight * (
        tables["hhe_cie"][temperature_index + 1] - tables["hhe_cie"][temperature_index]
    )
    metal_cie = tables["metal_cie"][temperature_index] + temperature_weight * (
        tables["metal_cie"][temperature_index + 1] - tables["metal_cie"][temperature_index]
    )
    lambda_cie_table = metallicity * metal_cie + primordial_cie
    low_temperature_e1 = math.exp(-1.184e5 / (temperature + 1.0e3))
    low_temperature_e2 = math.exp(-92.0 / temperature)
    lambda_low_temperature = metallicity * (
        2.0e-19 * low_temperature_e1
        + 2.8e-28 * math.sqrt(temperature) * low_temperature_e2
    )
    lambda_cie = temperature_in_table * lambda_cie_table + (
        1.0 - temperature_in_table
    ) * low_temperature * lambda_low_temperature
    return lambda_cie, lambda_pie


def net_cooling_rate(
    x: float, y: float, z: float, dx: float, rho: float, temperature: float, scalar0: float
) -> float:
    n_h = rho * NH_PER_DENSITY
    metallicity = scalar0 / COOLING_SOLAR_METALLICITY
    lambda_cie, lambda_pie = cooling_lambdas(temperature, n_h, metallicity)
    radius_squared = x * x + y * y
    radius = math.sqrt(radius_squared)
    horizontal_falloff = math.exp(-radius / COOLING_HSCALE_RADIUS)
    vertical_scale_squared = COOLING_HSCALE_HEIGHT**2 * (
        1.0 + radius_squared / COOLING_HSCALE_RADIUS**2
    )
    vertical_falloff = math.exp(-(z * z) / vertical_scale_squared)
    gamma_heating = (
        COOLING_HRATE
        * COOLING_HSCALE_NORM
        * NH_PER_DENSITY
        * horizontal_falloff
        * vertical_falloff
    )
    hot = 1.0 if temperature > 1.0e4 else 0.0
    inverse_ratio = 1.0e4 / temperature
    damping_factor = hot * inverse_ratio**8 + (1.0 - hot)
    gamma_heating *= damping_factor
    neutral_fraction = 1.0 - 0.5 * (1.0 + math.tanh((temperature - 8.0e3) / 1.5e3))
    tau = neutral_fraction * n_h * 1.0e-17 * dx * LENGTH_CGS
    pie_fraction = math.exp(-tau)
    lambda_cooling = (1.0 - pie_fraction) * lambda_cie + pie_fraction * lambda_pie
    return n_h * (n_h * lambda_cooling - gamma_heating)


def derived_values(
    block: BlockSpec, cell: int, stored: Sequence[float], gamma: float = GAMMA_DEFAULT
) -> Tuple[Mapping[str, float], Mapping[str, float]]:
    rho, vx, vy, vz, eint, scalar0, scalar1, scalar2 = stored
    x, y, z, volume = cell_center(block, cell)
    cylindrical_radius = math.sqrt(x * x + y * y)
    radius = math.sqrt(cylindrical_radius * cylindrical_radius + z * z)
    costheta = z / radius if radius > 0.0 else 1.0
    velocity_r = (vx * x + vy * y + vz * z) / radius if radius > 0.0 else 0.0
    velocity_theta = (
        z * (vx * x + vy * y) / (radius * cylindrical_radius)
        - vz * cylindrical_radius / radius
        if radius > 0.0 and cylindrical_radius > 0.0
        else 0.0
    )
    velocity_phi = (-vx * y + vy * x) / cylindrical_radius if cylindrical_radius > 0.0 else 0.0
    velocity_squared = vx * vx + vy * vy + vz * vz
    speed = math.sqrt(velocity_squared)
    temperature = eint / rho
    sound_speed = math.sqrt(gamma * (gamma - 1.0) * eint / rho)
    edot_kin = 0.5 * rho * velocity_squared * velocity_r
    edot_th = gamma * eint * velocity_r
    dust_fraction = scalar1 + scalar2
    total_metal_fraction = scalar0 + dust_fraction
    vertical_velocity_outward = vz if z >= 0.0 else -vz
    enthalpy_plus_ke = 0.5 * rho * velocity_squared + gamma * eint
    temperature_kelvin = temperature * TEMPERATURE_K_PER_RATIO
    cooling_rate = net_cooling_rate(
        x,
        y,
        z,
        (block.geometry[1] - block.geometry[0]) / block.shape[0],
        rho,
        temperature_kelvin,
        scalar0,
    )
    cooling_time_myr = eint * PRESSURE_CGS / (
        cooling_rate + FLOAT32_TINY * PRESSURE_CGS / TIME_CGS
    ) / MYR_CGS
    cooling_luminosity_erg_s = cooling_rate * volume * LENGTH_CGS**3
    values = {
        "coord_r": radius,
        "coord_costheta": costheta,
        "coord_abscostheta": abs(costheta),
        "hydro_w_d": rho,
        "temperature": temperature,
        "hydro_w_e": eint,
        "vel_sph_r": velocity_r,
        "vel_sph_theta": velocity_theta,
        "vel_sph_phi": velocity_phi,
        "hydro_w_s_00": scalar0,
        "hydro_w_s_01": scalar1,
        "hydro_w_s_02": scalar2,
        "speed": speed,
        "mach": speed / sound_speed,
        "radial_mach": velocity_r / sound_speed,
        "entropy_proxy": (gamma - 1.0) * eint / rho**gamma,
        "specific_angular_momentum_z": x * vy - y * vx,
        "specific_angular_momentum_z_kpc_km_s": (x * vy - y * vx)
        * VELOCITY_KM_S,
        "temperature_kelvin": temperature_kelvin,
        "hydrogen_number_density_cm3": rho * NH_PER_DENSITY,
        "pressure_over_kb_k_cm3": eint * PRESSURE_OVER_KB_PER_EINT,
        "velocity_r_km_s": velocity_r * VELOCITY_KM_S,
        "velocity_theta_km_s": velocity_theta * VELOCITY_KM_S,
        "velocity_phi_km_s": velocity_phi * VELOCITY_KM_S,
        "absolute_velocity_r_km_s": abs(velocity_r) * VELOCITY_KM_S,
        "sound_speed_km_s": sound_speed * VELOCITY_KM_S,
        "absolute_radial_mach": abs(velocity_r / sound_speed),
        "absolute_z": abs(z),
        "cylindrical_radius": cylindrical_radius,
        "dust_to_total_metal": (
            dust_fraction / total_metal_fraction
            if total_metal_fraction > 0.0
            else 0.0
        ),
        "small_grain_fraction": (
            scalar1 / dust_fraction if dust_fraction > 0.0 else 0.0
        ),
        "cooling_rate_erg_s_cm3": cooling_rate,
        "cooling_time_myr": cooling_time_myr,
    }
    weights = {
        "volume": volume,
        "mass": volume * rho,
        "mdot_sph": volume * rho * velocity_r,
        "mdot_sph_out": volume * rho * max(velocity_r, 0.0),
        "mdot_sph_in": volume * rho * min(velocity_r, 0.0),
        "edot_sph": volume * (edot_kin + edot_th),
        "edot_sph_out": volume * ((edot_kin + edot_th) if velocity_r > 0.0 else 0.0),
        "edot_sph_in": volume * ((edot_kin + edot_th) if velocity_r < 0.0 else 0.0),
        "edot_sph_kin": volume * edot_kin,
        "edot_sph_th": volume * edot_th,
        "mdot_sph_in_abs": volume * rho * max(-velocity_r, 0.0),
        "edot_sph_in_abs": volume
        * (-(edot_kin + edot_th) if velocity_r < 0.0 else 0.0),
        "edot_sph_kin_out": volume * (edot_kin if velocity_r > 0.0 else 0.0),
        "edot_sph_kin_in_abs": volume
        * (-edot_kin if velocity_r < 0.0 else 0.0),
        "edot_sph_th_out": volume * (edot_th if velocity_r > 0.0 else 0.0),
        "edot_sph_th_in_abs": volume
        * (-edot_th if velocity_r < 0.0 else 0.0),
        "edot_cool": cooling_luminosity_erg_s,
        "radial_ram_out": volume * rho * max(velocity_r, 0.0) ** 2,
        "vertical_mdot_out": volume
        * rho
        * max(vertical_velocity_outward, 0.0),
        "vertical_mdot_in_abs": volume
        * rho
        * max(-vertical_velocity_outward, 0.0),
        "vertical_edot_out": volume
        * enthalpy_plus_ke
        * max(vertical_velocity_outward, 0.0),
        "vertical_ram_out": volume
        * rho
        * max(vertical_velocity_outward, 0.0) ** 2,
        "gas_metal_mdot_out": volume
        * rho
        * max(velocity_r, 0.0)
        * scalar0,
        "gas_metal_mdot_in_abs": volume
        * rho
        * max(-velocity_r, 0.0)
        * scalar0,
        "total_metal_mdot_out": volume
        * rho
        * max(velocity_r, 0.0)
        * total_metal_fraction,
        "total_metal_mdot_in_abs": volume
        * rho
        * max(-velocity_r, 0.0)
        * total_metal_fraction,
        "dust_mdot_out": volume
        * rho
        * max(velocity_r, 0.0)
        * dust_fraction,
        "dust_mdot_in_abs": volume
        * rho
        * max(-velocity_r, 0.0)
        * dust_fraction,
        "total_metal_mass": volume * rho * total_metal_fraction,
        "dust_mass": volume * rho * dust_fraction,
    }
    return values, weights


def reference_histograms(
    blocks: Sequence[BlockSpec], products: Iterable[Product], gamma: float = GAMMA_DEFAULT
) -> Dict[str, Dict[int, float]]:
    """Return sparse flat-index histograms without using the reducer's headers."""

    product_list = tuple(products)
    result: Dict[str, Dict[int, float]] = {product.identifier: {} for product in product_list}
    for block, cell, stored in iter_cells(blocks):
        values, weights = derived_values(block, cell, stored, gamma)
        for product in product_list:
            flat = sum(
                bin_index(values[axis.variable], axis) * stride
                for axis, stride in zip(product.axes, product.strides)
            )
            histogram = result[product.identifier]
            histogram[flat] = histogram.get(flat, 0.0) + weights[product.weight]
    return result


def parse_header(path: Path) -> Dict[str, object]:
    """Small independent parser used to validate the reducer's emitted metadata."""

    raw = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#") and "=" in stripped:
            key, value = stripped.split("=", 1)
            raw[key.strip()] = value.strip()
    return raw
