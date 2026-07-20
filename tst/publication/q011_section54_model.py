"""Pure, bounded Q-011 Section 5.4 model helpers.

This module freezes locally decidable model assumptions only.  It does not
authorize a launch, inspect campaign policy, or qualify Section 5.4 evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real
import re
from typing import Mapping


INTEGRATOR = "vl2"
RUNTIME_MOMENTUM_STATE = "momentum_p_over_m"
DECK_INITIAL_STATE = "momentum"
LIGHT_SPEED = 10000.0
UPSTREAM_SPEED_U0 = 30.0
IDEAL_SURFACE_SPEED = 10.0
RESTART_SCHEMA = 8
BASE_DECK_PATH = (
    "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
)
FLOAT32_PROJECTION_UNCERTAINTY = (
    "Particle-VTK physical velocity components are float32 projections. "
    "Reconstructed chi is therefore a float32-projection estimate, not an "
    "exact recovery of pre-serialization momentum; uncertainty is amplified "
    "as velocity approaches C."
)

_RUNTIME_MODEL_PREFIX = "PIC runtime model:"
_RUNTIME_TOKEN_PATTERN = re.compile(r"([A-Za-z][A-Za-z0-9_]*)=([^\s=]+)")
_INTEGER_PATTERN = re.compile(r"(?:0|[1-9][0-9]*)")
_OVERRIDE_PATTERN = re.compile(r"([^/=\s]+)/([^/=\s]+)=([^\s=]+)")


class ModelContractError(ValueError):
    """Raised when a bounded Section 5.4 model assumption is not satisfied."""


@dataclass(frozen=True)
class VariantBinding:
    """Canonical base-deck and model-override binding for one grid variant."""

    variant: str
    deck_path: str
    model_launch_overrides: tuple[str, ...]


@dataclass(frozen=True)
class RuntimeIdentity:
    """Typed projection of the required PIC runtime-model identity."""

    integrator: str
    state: str
    light_speed: float
    restart_schema: int
    fields: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class DeckContract:
    """Typed projection of the deck assumptions used by these helpers."""

    variant: str
    integrator: str
    initial_state: str
    light_speed: float
    gamma: float
    upstream_speed_u0: float
    shock_speed_model: str
    nx1: int
    nx2: int
    nx3: int
    refinement: str
    num_levels: int
    curvature_amr: bool
    root_dx: float
    root_dy: float
    finest_dx: float
    finest_dy: float


_VARIANT_BINDINGS = {
    "coarse_uniform_dx12": VariantBinding(
        variant="coarse_uniform_dx12",
        deck_path=BASE_DECK_PATH,
        model_launch_overrides=(
            "mesh_refinement/refinement=none",
            "mesh_refinement/num_levels=1",
            "problem/ps_enable_curvature_amr=false",
        ),
    ),
    "three_level_amr_root_dx12_finest_dx3": VariantBinding(
        variant="three_level_amr_root_dx12_finest_dx3",
        deck_path=BASE_DECK_PATH,
        model_launch_overrides=(),
    ),
    "fine_uniform_dx3": VariantBinding(
        variant="fine_uniform_dx3",
        deck_path=BASE_DECK_PATH,
        model_launch_overrides=(
            "mesh/nx1=16000",
            "mesh/nx2=1040",
            "mesh_refinement/refinement=none",
            "mesh_refinement/num_levels=1",
            "problem/ps_enable_curvature_amr=false",
        ),
    ),
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ModelContractError(message)


def _finite_real(value: object, label: str) -> float:
    _require(
        isinstance(value, Real) and not isinstance(value, bool),
        f"{label}: expected a finite real number",
    )
    try:
        measured = float(value)
    except (OverflowError, TypeError, ValueError) as error:
        raise ModelContractError(f"{label}: expected a finite real number") from error
    _require(math.isfinite(measured), f"{label}: expected a finite real number")
    return measured


def _finite_float_text(value: str, label: str) -> float:
    try:
        measured = float(value)
    except (TypeError, ValueError) as error:
        raise ModelContractError(f"{label}: expected a finite number") from error
    _require(math.isfinite(measured), f"{label}: expected a finite number")
    return measured


def _positive_int_text(value: str, label: str) -> int:
    _require(
        _INTEGER_PATTERN.fullmatch(value) is not None,
        f"{label}: expected a positive integer",
    )
    measured = int(value)
    _require(measured > 0, f"{label}: expected a positive integer")
    return measured


def _bool_text(value: str, label: str) -> bool:
    lowered = value.lower()
    if lowered in {"true", "1"}:
        return True
    if lowered in {"false", "0"}:
        return False
    raise ModelContractError(f"{label}: expected a boolean")


def _require_close(measured: float, expected: float, label: str) -> None:
    _require(
        math.isclose(measured, expected, rel_tol=1.0e-11, abs_tol=1.0e-12),
        f"{label}: expected {expected!r}, measured {measured!r}",
    )


def canonical_variant_identity(variant: object) -> str:
    """Return an exact bounded variant identity, rejecting aliases."""
    _require(type(variant) is str, "variant: expected a canonical string identity")
    _require(
        variant in _VARIANT_BINDINGS,
        f"variant: unsupported canonical identity {variant!r}",
    )
    return variant


def variant_binding(variant: object) -> VariantBinding:
    """Return the immutable base-deck and model-override binding for a variant."""
    return _VARIANT_BINDINGS[canonical_variant_identity(variant)]


def _parse_model_launch_overrides(
    overrides: object,
) -> tuple[tuple[str, ...], tuple[tuple[str, str, str], ...]]:
    _require(
        not isinstance(overrides, (str, bytes)),
        "model launch overrides: expected an iterable of override records",
    )
    try:
        records = tuple(overrides)  # type: ignore[arg-type]
    except TypeError as error:
        raise ModelContractError(
            "model launch overrides: expected an iterable of override records"
        ) from error
    parsed = []
    names = set()
    for record in records:
        _require(type(record) is str, "model launch override: expected a string")
        match = _OVERRIDE_PATTERN.fullmatch(record)
        _require(match is not None, f"model launch override: malformed record {record!r}")
        block, name, value = match.groups()
        key = (block, name)
        _require(key not in names, f"model launch override: duplicate {block}/{name}")
        names.add(key)
        parsed.append((block, name, value))
    return records, tuple(parsed)


def require_variant_binding(
    variant: object,
    deck_path: object,
    model_launch_overrides: object,
) -> VariantBinding:
    """Require the exact base deck and model overrides for one canonical variant."""
    binding = variant_binding(variant)
    _require(type(deck_path) is str, "deck path: expected a string")
    _require(
        deck_path == binding.deck_path,
        f"{binding.variant}: deck path differs from the frozen binding",
    )
    records, _ = _parse_model_launch_overrides(model_launch_overrides)
    _require(
        records == binding.model_launch_overrides,
        f"{binding.variant}: model launch overrides differ from the frozen binding",
    )
    return binding


def _parse_athinput(deck_text: object) -> dict[str, dict[str, str]]:
    _require(type(deck_text) is str, "deck: expected UTF-8 text")
    blocks: dict[str, dict[str, str]] = {}
    current = None
    for lineno, raw_line in enumerate(deck_text.splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            _require(line.endswith(">"), f"deck:{lineno}: malformed block header")
            block = line[1:-1].strip()
            _require(bool(block), f"deck:{lineno}: empty block name")
            _require(block not in blocks, f"deck:{lineno}: duplicate block {block}")
            blocks[block] = {}
            current = block
            continue
        _require(current is not None, f"deck:{lineno}: parameter before first block")
        _require("=" in line, f"deck:{lineno}: malformed parameter line")
        name, value = (item.strip() for item in line.split("=", 1))
        _require(bool(name) and bool(value), f"deck:{lineno}: empty parameter or value")
        _require(
            name not in blocks[current],
            f"deck:{lineno}: duplicate {current}/{name}",
        )
        blocks[current][name] = value
    return blocks


def _parameter(blocks: Mapping[str, Mapping[str, str]], block: str, name: str) -> str:
    try:
        return blocks[block][name]
    except KeyError as error:
        raise ModelContractError(f"deck: missing {block}/{name}") from error


def _apply_model_launch_overrides(
    blocks: Mapping[str, Mapping[str, str]],
    overrides: object,
) -> dict[str, dict[str, str]]:
    _, parsed = _parse_model_launch_overrides(overrides)
    materialized = {block: dict(parameters) for block, parameters in blocks.items()}
    for block, name, value in parsed:
        _require(block in materialized, f"model launch override: missing block {block}")
        _require(
            name in materialized[block],
            f"model launch override: missing parameter {block}/{name}",
        )
        materialized[block][name] = value
    return materialized


def _typed_deck_contract(blocks: Mapping[str, Mapping[str, str]]) -> DeckContract:
    integrator = _parameter(blocks, "time", "integrator")
    _require(
        integrator == INTEGRATOR,
        "deck time/integrator: expected vl2",
    )
    for name, expected in (
        ("pusher", "boris_tsc"),
        ("deposit_moments", "true"),
        ("deposit_order", "2"),
        ("couple_moments_to_mhd", "true"),
        ("couple_moments_momentum_to_mhd", "true"),
        ("couple_moments_energy_to_mhd", "true"),
        ("pic_background_mode", "coupled"),
        ("pic_feedback_mode", "coupled"),
    ):
        _require(
            _parameter(blocks, "particles", name) == expected,
            f"deck particles/{name}: expected {expected}",
        )
    initial_state = _parameter(blocks, "particles", "pic_cr_initial_state")
    _require(
        initial_state == DECK_INITIAL_STATE,
        "deck particles/pic_cr_initial_state: expected momentum",
    )
    light_speed = _finite_float_text(
        _parameter(blocks, "particles", "pic_cr_light_speed"),
        "deck particles/pic_cr_light_speed",
    )
    _require_close(light_speed, LIGHT_SPEED, "deck particles/pic_cr_light_speed")

    gamma = _finite_float_text(_parameter(blocks, "mhd", "gamma"), "deck mhd/gamma")
    _require_close(gamma, 5.0 / 3.0, "deck mhd/gamma")
    upstream_speed = _finite_float_text(
        _parameter(blocks, "problem", "ps_u0"), "deck problem/ps_u0"
    )
    _require_close(upstream_speed, UPSTREAM_SPEED_U0, "deck problem/ps_u0")
    _require(
        _parameter(blocks, "problem", "pgen_name") == "pic_parallel_shock",
        "deck problem/pgen_name: expected pic_parallel_shock",
    )
    shock_speed_model = _parameter(blocks, "problem", "ps_shock_speed_model")
    _require(
        shock_speed_model == "ideal_surface",
        "deck problem/ps_shock_speed_model: expected ideal_surface",
    )
    derived_surface_speed = 0.5 * (gamma - 1.0) * upstream_speed
    _require_close(derived_surface_speed, IDEAL_SURFACE_SPEED, "deck ideal surface speed")

    nx1 = _positive_int_text(_parameter(blocks, "mesh", "nx1"), "deck mesh/nx1")
    nx2 = _positive_int_text(_parameter(blocks, "mesh", "nx2"), "deck mesh/nx2")
    nx3 = _positive_int_text(_parameter(blocks, "mesh", "nx3"), "deck mesh/nx3")
    _require(nx3 == 1, "deck mesh/nx3: expected collapsed x3 dimension")
    x1min = _finite_float_text(_parameter(blocks, "mesh", "x1min"), "deck mesh/x1min")
    x1max = _finite_float_text(_parameter(blocks, "mesh", "x1max"), "deck mesh/x1max")
    x2min = _finite_float_text(_parameter(blocks, "mesh", "x2min"), "deck mesh/x2min")
    x2max = _finite_float_text(_parameter(blocks, "mesh", "x2max"), "deck mesh/x2max")
    x3min = _finite_float_text(_parameter(blocks, "mesh", "x3min"), "deck mesh/x3min")
    x3max = _finite_float_text(_parameter(blocks, "mesh", "x3max"), "deck mesh/x3max")
    for measured, expected, label in (
        (x1min, 0.0, "deck mesh/x1min"),
        (x1max, 48000.0, "deck mesh/x1max"),
        (x2min, 0.0, "deck mesh/x2min"),
        (x2max, 3120.0, "deck mesh/x2max"),
        (x3min, 0.0, "deck mesh/x3min"),
        (x3max, 1.0, "deck mesh/x3max"),
    ):
        _require_close(measured, expected, label)
    for name, expected in (("nx1", 20), ("nx2", 20), ("nx3", 1)):
        measured = _positive_int_text(
            _parameter(blocks, "meshblock", name), f"deck meshblock/{name}"
        )
        _require(measured == expected, f"deck meshblock/{name}: expected {expected}")

    refinement = _parameter(blocks, "mesh_refinement", "refinement")
    num_levels = _positive_int_text(
        _parameter(blocks, "mesh_refinement", "num_levels"),
        "deck mesh_refinement/num_levels",
    )
    curvature_amr = _bool_text(
        _parameter(blocks, "problem", "ps_enable_curvature_amr"),
        "deck problem/ps_enable_curvature_amr",
    )
    layout = (nx1, nx2, refinement, num_levels, curvature_amr)
    if layout == (4000, 260, "none", 1, False):
        variant = "coarse_uniform_dx12"
    elif layout == (4000, 260, "adaptive", 3, True):
        variant = "three_level_amr_root_dx12_finest_dx3"
    elif layout == (16000, 1040, "none", 1, False):
        variant = "fine_uniform_dx3"
    else:
        raise ModelContractError(f"deck: unsupported bounded grid layout {layout!r}")

    root_dx = (x1max - x1min) / nx1
    root_dy = (x2max - x2min) / nx2
    refinement_scale = 2 ** (num_levels - 1) if refinement == "adaptive" else 1
    finest_dx = root_dx / refinement_scale
    finest_dy = root_dy / refinement_scale
    expected_finest_spacing = 12.0 if variant == "coarse_uniform_dx12" else 3.0
    _require_close(finest_dx, expected_finest_spacing, "deck finest dx")
    _require_close(finest_dy, expected_finest_spacing, "deck finest dy")
    return DeckContract(
        variant=variant,
        integrator=integrator,
        initial_state=initial_state,
        light_speed=light_speed,
        gamma=gamma,
        upstream_speed_u0=upstream_speed,
        shock_speed_model=shock_speed_model,
        nx1=nx1,
        nx2=nx2,
        nx3=nx3,
        refinement=refinement,
        num_levels=num_levels,
        curvature_amr=curvature_amr,
        root_dx=root_dx,
        root_dy=root_dy,
        finest_dx=finest_dx,
        finest_dy=finest_dy,
    )


def parse_deck_contract(deck_text: object) -> DeckContract:
    """Parse a typed, bounded coarse-uniform, AMR, or fine-uniform deck contract."""
    return _typed_deck_contract(_parse_athinput(deck_text))


def parse_bound_variant_deck_contract(
    variant: object,
    deck_path: object,
    deck_text: object,
    model_launch_overrides: object,
) -> DeckContract:
    """Materialize and validate one exact variant binding from the AMR base deck."""
    binding = require_variant_binding(variant, deck_path, model_launch_overrides)
    blocks = _parse_athinput(deck_text)
    base_contract = _typed_deck_contract(blocks)
    _require(
        base_contract.variant == "three_level_amr_root_dx12_finest_dx3",
        "base deck: expected the frozen AMR layout",
    )
    contract = _typed_deck_contract(
        _apply_model_launch_overrides(blocks, binding.model_launch_overrides)
    )
    _require(
        contract.variant == binding.variant,
        f"{binding.variant}: materialized deck variant does not match its binding",
    )
    return contract


def parse_runtime_identity_line(line: object) -> RuntimeIdentity:
    """Parse one PIC runtime-model line and require the Section 5.4 identity."""
    _require(type(line) is str, "runtime identity: expected text")
    _require("\n" not in line and "\r" not in line, "runtime identity: expected one line")
    prefix = _RUNTIME_MODEL_PREFIX + " "
    _require(line.startswith(prefix), "runtime identity: malformed prefix")
    tokens = line[len(prefix):].split(" ")
    _require(bool(tokens) and all(tokens), "runtime identity: malformed token spacing")
    fields = {}
    for token in tokens:
        match = _RUNTIME_TOKEN_PATTERN.fullmatch(token)
        _require(match is not None, f"runtime identity: malformed token {token!r}")
        name, value = match.groups()
        _require(name not in fields, f"runtime identity: duplicate field {name}")
        fields[name] = value
    required = {
        "integrator",
        "state",
        "C",
        "background",
        "feedback",
        "induction",
        "deposition",
        "restart_schema",
    }
    missing = sorted(required - set(fields))
    _require(not missing, f"runtime identity: missing fields {', '.join(missing)}")
    _require(
        fields["integrator"] == INTEGRATOR,
        "runtime identity integrator: expected vl2",
    )
    _require(
        fields["state"] == RUNTIME_MOMENTUM_STATE,
        "runtime identity state: expected momentum_p_over_m",
    )
    for name, expected in (
        ("background", "coupled"),
        ("feedback", "coupled"),
        ("induction", "ideal_mhd_only"),
        ("deposition", "tsc"),
    ):
        _require(
            fields[name] == expected,
            f"runtime identity {name}: expected {expected}",
        )
    light_speed = _finite_float_text(fields["C"], "runtime identity C")
    _require_close(light_speed, LIGHT_SPEED, "runtime identity C")
    _require(
        _INTEGER_PATTERN.fullmatch(fields["restart_schema"]) is not None,
        "runtime identity restart_schema: expected an integer",
    )
    restart_schema = int(fields["restart_schema"])
    _require(
        restart_schema == RESTART_SCHEMA,
        f"runtime identity restart_schema: expected {RESTART_SCHEMA}",
    )
    return RuntimeIdentity(
        integrator=fields["integrator"],
        state=fields["state"],
        light_speed=light_speed,
        restart_schema=restart_schema,
        fields=tuple(fields.items()),
    )


def parse_runtime_identity(stdout: object) -> RuntimeIdentity:
    """Require exactly one Section 5.4 runtime-model line in Athena stdout."""
    _require(type(stdout) is str, "runtime stdout: expected text")
    lines = [
        line for line in stdout.splitlines() if line.startswith(_RUNTIME_MODEL_PREFIX)
    ]
    _require(
        len(lines) == 1,
        "runtime stdout: expected exactly one PIC runtime-model identity line",
    )
    return parse_runtime_identity_line(lines[0])


def x_ideal(time: object) -> float:
    """Return the frozen Section 5.4 ideal shock position, x_ideal(t) = 10*t."""
    measured_time = _finite_real(time, "time")
    _require(measured_time >= 0.0, "time: expected a nonnegative value")
    position = IDEAL_SURFACE_SPEED * measured_time
    _require(math.isfinite(position), "ideal shock position: expected a finite value")
    return position


def reconstruct_chi_from_physical_speed(
    physical_speed: object,
) -> float:
    """Reconstruct chi from a physical speed magnitude under the frozen model.

    Particle-VTK velocities are serialized float32 projections.  The returned
    value is consequently a projection estimate, not exact pre-serialization
    momentum recovery.
    """
    speed = _finite_real(physical_speed, "physical speed")
    _require(speed >= 0.0, "physical speed: expected a nonnegative magnitude")
    _require(speed < LIGHT_SPEED, "physical speed: expected a subluminal magnitude")
    speed_squared = speed * speed
    denominator = 1.0 - speed_squared / (LIGHT_SPEED * LIGHT_SPEED)
    _require(
        math.isfinite(denominator) and denominator > 0.0,
        "chi reconstruction: invalid relativistic denominator",
    )
    chi = (speed_squared / denominator) / (UPSTREAM_SPEED_U0 * UPSTREAM_SPEED_U0)
    _require(math.isfinite(chi), "chi reconstruction: expected a finite value")
    return chi


def reconstruct_chi_from_physical_velocity(
    velocity_xyz: object,
) -> float:
    """Reconstruct chi from one float32-projected physical velocity vector."""
    _require(
        not isinstance(velocity_xyz, (str, bytes)),
        "physical velocity: expected exactly three components",
    )
    try:
        components = tuple(velocity_xyz)  # type: ignore[arg-type]
    except TypeError as error:
        raise ModelContractError(
            "physical velocity: expected exactly three components"
        ) from error
    _require(len(components) == 3, "physical velocity: expected exactly three components")
    vx, vy, vz = (
        _finite_real(component, f"physical velocity component {index}")
        for index, component in enumerate(components)
    )
    speed = math.hypot(vx, vy, vz)
    _require(math.isfinite(speed), "physical velocity: expected a finite magnitude")
    return reconstruct_chi_from_physical_speed(speed)
