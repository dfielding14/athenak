from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
BASE_DECK = REPO_ROOT / "inputs/tests/pic_static_load_balance_restart.athinput"
OVERRIDE_DECK = (
    REPO_ROOT / "inputs/tests/pic_static_load_balance_restart_override.athinput"
)
STATIC_INTERVAL = "pic_static_load_balance_interval"
STATIC_COST = "pic_static_load_balance_cost_per_particle"
PHYSICAL_COST = "pic_load_balance_cost_per_particle"


def _function_body(source: str, signature: str) -> str:
    start = source.index(signature)
    brace = source.index("{", start)
    depth = 0
    for index in range(brace, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[brace:index + 1]
    raise AssertionError(f"unterminated function: {signature}")


def _parse_deck(path: Path) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    block = ""
    for raw_line in path.read_text(encoding="ascii").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1]
            blocks.setdefault(block, {})
            continue
        name, value = line.split("=", 1)
        blocks[block][name.strip()] = value.strip()
    return blocks


def test_static_cadence_is_disabled_by_default_and_gated_to_nonadaptive_pic() -> None:
    driver = (REPO_ROOT / "src/driver/driver.cpp").read_text(encoding="ascii")
    gate = _function_body(driver, "bool Driver::StaticParticleLoadBalanceDue")
    assert "pic_static_load_balance_interval_(0)" in driver
    assert "pic_static_load_balance_cost_per_particle_(0.0)" in driver
    assert "pic_static_load_balance_interval_ > 0" in gate
    assert "pic_static_load_balance_cost_per_particle_ > 0.0" in gate
    assert "!pm->adaptive" in gate
    assert "pm->pmb_pack->ppart != nullptr" in gate
    assert "pm->ncycle % pic_static_load_balance_interval_ == 0" in gate

    execute = _function_body(driver, "void Driver::Execute(")
    assert "} else if (StaticParticleLoadBalanceDue(pmesh))" in execute
    assert "RedistributeStaticMeshBlocks" in execute


def test_static_weight_is_not_particle_restart_physical_metadata() -> None:
    particles = (REPO_ROOT / "src/particles/particles.hpp").read_text(encoding="ascii")
    metadata = _function_body(particles, "void FillRestartModelMetadata")
    assert PHYSICAL_COST in metadata
    assert STATIC_COST not in metadata
    assert STATIC_INTERVAL not in particles
    assert STATIC_COST not in particles

    refinement = (
        REPO_ROOT / "src/mesh/mesh_refinement.cpp"
    ).read_text(encoding="ascii")
    adaptive = _function_body(
        refinement,
        "void MeshRefinement::RedistAndRefineMeshBlocks(ParameterInput *pin, "
        "int nnew, int ndel)",
    )
    static = _function_body(
        refinement, "void MeshRefinement::RedistributeStaticMeshBlocks"
    )
    assert f"ppart->{PHYSICAL_COST}" in adaptive
    assert "pdriver, pin, 0, 0, cost_per_particle" in static

    reinitialize = _function_body(
        refinement, "void MeshRefinement::ReinitializeAfterRedistribution"
    )
    assert "InitBoundaryValuesAndPrimitives" in reinitialize
    for module in ("phydro", "pmhd", "prad", "pz4c"):
        assert f"pmbp->{module}->NewTimeStep" in reinitialize
    assert "AssignParticleAwareCosts(new_cost_eachmb, new_nmb, cost_per_particle)" in (
        refinement
    )


def test_static_controls_validate_and_reserve_capacity_only_when_enabled() -> None:
    driver = (REPO_ROOT / "src/driver/driver.cpp").read_text(encoding="ascii")
    assert f'"<particles>/{STATIC_INTERVAL} must be >= 0"' in driver
    assert f'"<particles>/{STATIC_COST} must be "' in driver
    assert "finite and >= 0" in driver

    build_tree = (REPO_ROOT / "src/mesh/build_tree.cpp").read_text(encoding="ascii")
    request = _function_body(build_tree, "bool StaticParticleLoadBalanceRequested")
    assert STATIC_INTERVAL in request
    assert STATIC_COST in request
    assert "!adaptive && StaticParticleLoadBalanceRequested(pin)" in build_tree
    assert "if (multilevel || static_particle_lb)" in build_tree
    assert "static particle load balancing" in build_tree
    assert "max_nmb_per_rank" in build_tree


def test_uniform_shock_fast_path_preserves_exact_ledger_cost_guard() -> None:
    source = (
        REPO_ROOT / "src/pgen/tests/pic_parallel_shock.cpp"
    ).read_text(encoding="ascii")
    topology = _function_body(
        source, "bool ParallelShockMeshStateIsFixedUniform"
    )
    exact = _function_body(
        source, "bool ParallelShockExactMeshStateIsFixedUniform"
    )
    preparation = _function_body(
        source, "void PrepareParallelShockInjectionTransaction"
    )
    assert "pmesh->cost_eachmb[gid] <= 0.0F" in topology
    assert "require_unit_costs && pmesh->cost_eachmb[gid] != 1.0F" in topology
    assert "ParallelShockMeshStateIsFixedUniform(pmesh, true)" in exact
    assert "ParallelShockMeshStateIsFixedUniform(pm, false)" in preparation
    assert "static_particle_redistribution_enabled" in source
    assert '"pic_static_load_balance_interval", 0) > 0' in source
    assert '"pic_static_load_balance_cost_per_particle", 0.0) > 0.0' in source
    assert "!static_particle_redistribution_enabled" in source


def test_restart_override_keeps_physical_cost_zero() -> None:
    base = _parse_deck(BASE_DECK)
    override = _parse_deck(OVERRIDE_DECK)
    assert base["mesh_refinement"]["refinement"] == "static"
    assert int(base["mesh_refinement"]["max_nmb_per_rank"]) > 0
    assert float(base["particles"][PHYSICAL_COST]) == 0.0
    assert STATIC_INTERVAL not in base["particles"]
    assert STATIC_COST not in base["particles"]
    assert PHYSICAL_COST not in override["particles"]
    assert int(override["particles"][STATIC_INTERVAL]) > 0
    assert float(override["particles"][STATIC_COST]) > 0.0


def _executable() -> Path:
    executable_dir = os.environ.get("ATHENA_STATIC_PIC_LB_EXE_DIR")
    if not executable_dir:
        pytest.skip("ATHENA_STATIC_PIC_LB_EXE_DIR is required for runtime checks")
    executable = Path(executable_dir) / "athena"
    if not executable.is_file():
        raise RuntimeError(f"Athena executable is missing: {executable}")
    return executable


def _run(executable: Path, arguments: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(executable), *arguments],
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def _telemetry(stdout: str, name: str) -> float:
    match = re.search(rf"^{re.escape(name)}=([^\s]+)$", stdout, re.MULTILINE)
    if match is None:
        raise AssertionError(f"missing telemetry scalar: {name}")
    return float(match.group(1))


def test_runtime_restart_accepts_performance_only_static_cost(tmp_path: Path) -> None:
    executable = _executable()
    initial_dir = tmp_path / "initial"
    initial_dir.mkdir()
    initial = _run(executable, ["-i", str(BASE_DECK), "-d", str(initial_dir)])
    assert initial.returncode == 0, initial.stdout

    restarts = sorted((initial_dir / "rst").glob("*.rst"))
    assert restarts
    restart = restarts[-1]
    assert Path(str(restart) + ".manifest").is_file()

    continuation_dir = tmp_path / "continuation"
    continuation_dir.mkdir()
    continuation = _run(
        executable,
        [
            "-r",
            str(restart),
            "-i",
            str(OVERRIDE_DECK),
            "-d",
            str(continuation_dir),
        ],
    )
    assert continuation.returncode == 0, continuation.stdout
    assert "Particle restart physical-model metadata mismatch" not in continuation.stdout
    assert "cycle=2" in continuation.stdout
    initial_nmb = _telemetry(initial.stdout, "q017.telemetry.meshblocks.total")
    continuation_nmb = _telemetry(
        continuation.stdout, "q017.telemetry.meshblocks.total"
    )
    continuation_cost = _telemetry(
        continuation.stdout, "q017.telemetry.load.cost.total"
    )
    assert continuation_nmb == initial_nmb
    assert continuation_cost > continuation_nmb


def test_runtime_uniform_mesh_accepts_opt_in_particle_redistribution(
    tmp_path: Path,
) -> None:
    executable = _executable()
    initial_dir = tmp_path / "uniform_initial"
    initial_dir.mkdir()
    initial = _run(
        executable,
        [
            "-i",
            str(BASE_DECK),
            "-d",
            str(initial_dir),
            "job/basename=pic_uniform_particle_lb",
            "mesh_refinement/refinement=none",
        ],
    )
    assert initial.returncode == 0, initial.stdout
    assert _telemetry(initial.stdout, "q017.telemetry.amr.enabled") == 0.0
    restarts = sorted((initial_dir / "rst").glob("*.rst"))
    assert restarts

    continuation_dir = tmp_path / "uniform_continuation"
    continuation_dir.mkdir()
    continuation = _run(
        executable,
        [
            "-r",
            str(restarts[-1]),
            "-i",
            str(OVERRIDE_DECK),
            "-d",
            str(continuation_dir),
        ],
    )
    assert continuation.returncode == 0, continuation.stdout
    assert _telemetry(continuation.stdout, "q017.telemetry.amr.enabled") == 0.0
    assert _telemetry(
        continuation.stdout, "q017.telemetry.load.cost.total"
    ) > _telemetry(
        continuation.stdout, "q017.telemetry.meshblocks.total"
    )


@pytest.mark.parametrize(
    ("parameter", "value", "message"),
    [
        (STATIC_INTERVAL, "-1", f"<particles>/{STATIC_INTERVAL} must be >= 0"),
        (STATIC_COST, "-1.0", f"<particles>/{STATIC_COST} must be"),
    ],
)
def test_runtime_rejects_negative_static_controls(
    tmp_path: Path, parameter: str, value: str, message: str
) -> None:
    executable = _executable()
    deck = tmp_path / f"negative_{parameter}.athinput"
    deck.write_text(
        BASE_DECK.read_text(encoding="ascii")
        + "\n<particles>\n"
        + f"{STATIC_INTERVAL} = 1\n"
        + f"{STATIC_COST} = 1.0\n"
        + f"{parameter} = {value}\n",
        encoding="ascii",
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    result = _run(executable, ["-i", str(deck), "-d", str(run_dir)])
    assert result.returncode != 0, result.stdout
    assert message in result.stdout
