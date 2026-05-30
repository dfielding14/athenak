"""CPU parser-contract tests for the fail-closed relativistic HC pusher."""

import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
RELATIVISTIC_INPUT = (
    ROOT / "inputs" / "particles" / "cr_tracer_relativistic_contract.athinput")
LEGACY_INPUT = ROOT / "inputs" / "particles" / "cr_tracer_boris_uniform.athinput"

VALID_RELATIVISTIC_FLAGS = []


def _render_input(tmp_path, input_file=RELATIVISTIC_INPUT, drop_parameters=(),
                  drop_blocks=(), inject_before_particles="",
                  inject_in_particles="", append_text=""):
    """Render a case input while keeping generated artifacts out of the repository."""
    contents = input_file.read_text()
    if drop_parameters:
        lines = contents.splitlines()
        contents = "\n".join(
            line for line in lines
            if line.split("=", 1)[0].strip() not in drop_parameters) + "\n"
    if drop_blocks:
        kept_lines = []
        dropping = False
        for line in contents.splitlines():
            stripped = line.strip()
            if stripped.startswith("<") and stripped.endswith(">"):
                dropping = stripped[1:-1] in drop_blocks
            if not dropping:
                kept_lines.append(line)
        contents = "\n".join(kept_lines) + "\n"
    contents = contents.replace(
        "\n<particles>\n", f"\n{inject_before_particles}<particles>\n", 1)
    contents = contents.replace(
        "\n<particles>\n", f"\n<particles>\n{inject_in_particles}", 1)
    contents += append_text
    case_input = tmp_path / "cr_relativistic_contract.athinput"
    case_input.write_text(contents)
    return case_input


def _run_capture(tmp_path, flags, input_file=RELATIVISTIC_INPUT,
                 drop_parameters=(), drop_blocks=(), inject_before_particles="",
                 inject_in_particles="", append_text=""):
    """Run the compiled CPU binary and retain fatal diagnostics for assertions."""
    input_file = _render_input(
        tmp_path, input_file, drop_parameters, drop_blocks, inject_before_particles,
        inject_in_particles, append_text)
    binary = Path("./athena").resolve()
    command = [str(binary), "-i", str(input_file)] + list(flags)
    return subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)


def _run_relativistic(tmp_path, extra_flags=(), drop_parameters=(),
                      drop_blocks=(), inject_before_particles="", inject_in_particles="",
                      append_text=""):
    return _run_capture(tmp_path, VALID_RELATIVISTIC_FLAGS + list(extra_flags),
                        drop_parameters=drop_parameters,
                        drop_blocks=drop_blocks,
                        inject_before_particles=inject_before_particles,
                        inject_in_particles=inject_in_particles,
                        append_text=append_text)


def _assert_fatal(result, substring):
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "### FATAL ERROR" in output
    assert substring.lower() in output.lower(), output


def _assert_fails_closed(result):
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "### FATAL ERROR" in output, output


@pytest.mark.parametrize("alpha_s", ["1.0", "-1.0"])
def test_relativistic_hc_valid_constructor_only_cases_cpu(tmp_path, alpha_s):
    """Signed nonzero alpha_s values are accepted when no push cycle is requested."""
    result = _run_relativistic(tmp_path, [f"particles/alpha_s={alpha_s}"])
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    ("missing_field", "expected_substring"),
    [
        ("c_model", "c_model"),
        ("alpha_s", "alpha_s"),
        ("relativistic_field_source", "relativistic_field_source"),
    ],
)
def test_relativistic_hc_requires_explicit_new_fields_cpu(
        tmp_path, missing_field, expected_substring):
    """The relativistic mode fails closed when a required new field is omitted."""
    result = _run_relativistic(tmp_path, drop_parameters=(missing_field,))
    _assert_fatal(result, expected_substring)


@pytest.mark.parametrize("c_model", ["0.0", "-1.0", "nan", "inf", "-inf"])
def test_relativistic_hc_rejects_invalid_c_model_cpu(tmp_path, c_model):
    """c_model must be finite and strictly positive."""
    result = _run_relativistic(tmp_path, [f"particles/c_model={c_model}"])
    _assert_fatal(result, "c_model")


@pytest.mark.parametrize("alpha_s", ["0.0", "nan", "inf", "-inf"])
def test_relativistic_hc_rejects_invalid_alpha_s_cpu(tmp_path, alpha_s):
    """alpha_s must be finite and nonzero, while either sign remains valid."""
    result = _run_relativistic(tmp_path, [f"particles/alpha_s={alpha_s}"])
    _assert_fatal(result, "alpha_s")


@pytest.mark.parametrize("interpolation", ["tsc", "lin", "lin_legacy"])
def test_relativistic_hc_requires_trilinear_interpolation_cpu(
        tmp_path, interpolation):
    """Legacy gathers remain available to Boris but are rejected by relativistic HC."""
    result = _run_relativistic(
        tmp_path, [f"particles/interpolation={interpolation}"])
    _assert_fatal(result, "interpolation")


def test_relativistic_hc_requires_explicit_interpolation_cpu(tmp_path):
    """The relativistic gather must be selected explicitly rather than defaulted."""
    result = _run_relativistic(tmp_path, drop_parameters=("interpolation",))
    _assert_fatal(result, "interpolation")


@pytest.mark.parametrize(
    ("extra_flags", "expected_substring"),
    [
        (["particles/nspecies=2"], "nspecies"),
        (["mesh/nx3=1", "meshblock/nx3=1"], "3d"),
        (["mesh/ix1_bc=outflow", "mesh/ox1_bc=outflow"], "periodic"),
        (["time/evolution=dynamic", "mhd/rsolver=hlld"], "kinematic"),
        (["time/evolution=static"], "kinematic"),
        (["problem/prtcl_rst_flag=1"], "legacy particle restart"),
    ],
    ids=["single-species-only", "three-dimensional-only", "periodic-only",
         "kinematic-only", "static-rejected", "legacy-restart-rejected"],
)
def test_relativistic_hc_rejects_unsupported_constructor_modes_cpu(
        tmp_path, extra_flags, expected_substring):
    """Unsupported Phase-2 constructor combinations abort with useful diagnostics."""
    result = _run_relativistic(tmp_path, extra_flags)
    _assert_fatal(result, expected_substring)


def test_relativistic_hc_rejects_unimplemented_field_source_cpu(tmp_path):
    """Only the prescribed Phase-2 test field is accepted."""
    result = _run_relativistic(
        tmp_path, ["particles/relativistic_field_source=mhd"])
    _assert_fatal(result, "relativistic_field_source")


def test_relativistic_hc_requires_mhd_block_cpu(tmp_path):
    """The prescribed-test constructor still requires explicit MHD plumbing."""
    result = _run_relativistic(tmp_path, drop_blocks=("mhd",))
    _assert_fatal(result, "requires an <mhd> block")


def test_relativistic_hc_rejects_stale_background_key_cpu(tmp_path):
    """The superseded frozen-background spelling must not be silently ignored."""
    result = _run_relativistic(
        tmp_path, inject_in_particles="relativistic_background = frozen\n")
    _assert_fatal(result, "relativistic_background is unsupported")


@pytest.mark.parametrize(
    ("extra_flags", "inject_before_particles", "append_text", "expected_substring"),
    [
        (["mhd/eos=isothermal"], "iso_sound_speed = 1.0\n", "", "eos = ideal"),
        ([], "", "\n<mhd_srcterms>\n", "coupled physics blocks"),
        ([], "ohmic_resistivity = 0.1\n", "", "mhd diffusion"),
        ([], "viscosity = 0.1\n", "", "mhd diffusion"),
        ([], "conductivity = 0.1\n", "", "mhd diffusion"),
        ([], "tdep_conductivity = 0.1\n", "", "mhd diffusion"),
        ([], "nscalars = 1\n", "", "passive scalars"),
        ([], "", "\n<turb_driving>\n", "coupled physics blocks"),
        ([], "", "\n<shearing_box>\nqshear = 1.5\nomega0 = 1.0\n",
         "coupled physics blocks"),
        (["time/evolution=dynamic", "mhd/rsolver=llf"], "",
         "\n<coord>\nspecial_rel = true\n", "Newtonian coordinates"),
        (["time/evolution=dynamic", "mhd/rsolver=llf"], "",
         "\n<coord>\ngeneral_rel = true\nminkowski = true\n", "Newtonian coordinates"),
    ],
    ids=["ideal-eos-only", "mhd-source-terms", "ohmic-resistivity", "viscosity",
         "conductivity", "time-dependent-conductivity", "passive-scalars",
         "active-turbulence-forcing", "shearing-box-and-orbital-advection",
         "sr-coordinates", "gr-coordinates"],
)
def test_relativistic_hc_rejects_unqualified_physics_cpu(
        tmp_path, extra_flags, inject_before_particles, append_text,
        expected_substring):
    """Every sampled-MHD widening remains fail-closed during parser qualification."""
    result = _run_relativistic(
        tmp_path, extra_flags, inject_before_particles=inject_before_particles,
        append_text=append_text)
    _assert_fatal(result, expected_substring)


@pytest.mark.parametrize(
    ("append_text", "case_name"),
    [
        ("\n<ion-neutral>\ndrag_coeff = 1.0\n", "ion-neutral"),
        ("\n<adm>\n", "dynamical-gr"),
    ],
)
def test_relativistic_hc_inherited_constructor_rejections_fail_closed_cpu(
        tmp_path, append_text, case_name):
    """Coupled modes rejected before particle construction must still fail closed."""
    del case_name
    result = _run_relativistic(tmp_path, append_text=append_text)
    _assert_fails_closed(result)


def test_relativistic_hc_positive_cycle_aborts_as_unimplemented_cpu(tmp_path):
    """Constructor success must not expose an accidental production push path."""
    result = _run_relativistic(tmp_path, ["time/nlim=1"])
    _assert_fatal(result, "not implemented")


def test_relativistic_hc_rejects_outputs_before_writing_files_cpu(tmp_path):
    """Parser-only cases must not serialize misleading legacy particle artifacts."""
    output = "\n<output1>\nfile_type = prst\ndt = 0.01\n"
    result = _run_relativistic(tmp_path, append_text=output)
    _assert_fatal(result, "until relativistic outputs are qualified")
    assert not (tmp_path / "prst").exists()


def test_legacy_pushers_reject_relativistic_only_keys_cpu(tmp_path):
    """Legacy pushers must not silently ignore staged relativistic parameters."""
    flags = [
        "particles/pusher=boris",
        "time/evolution=dynamic",
        "mhd/rsolver=hlld",
        "time/nlim=1",
    ]
    result = _run_capture(tmp_path, flags)
    _assert_fatal(result, "relativistic-only")


def test_legacy_pushers_reject_stale_relativistic_background_cpu(tmp_path):
    """Legacy modes must also reject the superseded relativistic spelling."""
    flags = [
        "particles/pusher=boris",
        "time/evolution=dynamic",
        "mhd/rsolver=hlld",
        "time/nlim=1",
    ]
    result = _run_capture(
        tmp_path, flags, inject_in_particles="relativistic_background = frozen\n")
    _assert_fatal(result, "relativistic-only")


@pytest.mark.parametrize(
    ("pusher", "interpolation"),
    [
        ("drift", "tsc"),
        ("boris", "tsc"),
        ("boris", "lin"),
        ("boris", "lin_legacy"),
        ("boris", "trilinear"),
    ],
)
def test_legacy_particle_constructor_controls_remain_valid_cpu(
        tmp_path, pusher, interpolation):
    """Adding relativistic HC must not tighten existing drift or Boris parsing."""
    flags = [
        f"particles/pusher={pusher}",
        f"particles/interpolation={interpolation}",
        "time/nlim=0",
    ]
    result = _run_capture(tmp_path, flags, input_file=LEGACY_INPUT)
    assert result.returncode == 0, result.stdout + result.stderr
