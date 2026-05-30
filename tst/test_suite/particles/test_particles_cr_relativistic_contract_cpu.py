"""CPU parser-contract tests for the fail-closed relativistic HC pusher."""

import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
RELATIVISTIC_INPUT = (
    ROOT / "inputs" / "particles" / "cr_tracer_relativistic_contract.athinput")
LEGACY_INPUT = ROOT / "inputs" / "particles" / "cr_tracer_boris_uniform.athinput"

VALID_RELATIVISTIC_FLAGS = []
COUPLED_RELATIVISTIC_FLAGS = [
    "particles/relativistic_field_source=mhd_ideal",
    "time/evolution=dynamic",
    "mhd/rsolver=hlld",
]
COUPLED_DROP_PARAMETERS = ("cE0x", "cE0y", "cE0z")


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


def _run_coupled(tmp_path, extra_flags=(), drop_parameters=(), drop_blocks=(),
                 inject_before_particles="", inject_in_particles="", append_text="",
                 temporal_mode="frozen_tn"):
    temporal_input = (
        f"relativistic_temporal_sampling = {temporal_mode}\n"
        if temporal_mode is not None else "")
    return _run_capture(
        tmp_path, COUPLED_RELATIVISTIC_FLAGS + list(extra_flags),
        drop_parameters=COUPLED_DROP_PARAMETERS + tuple(drop_parameters),
        drop_blocks=drop_blocks, inject_before_particles=inject_before_particles,
        inject_in_particles=temporal_input + inject_in_particles,
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
        ("relativistic_initial_state", "relativistic_initial_state"),
        ("B0x", "B0x"),
        ("B0y", "B0y"),
        ("B0z", "B0z"),
        ("cE0x", "cE0x"),
        ("cE0y", "cE0y"),
        ("cE0z", "cE0z"),
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
    """Unqualified and CT-edge field-source spellings remain closed."""
    result = _run_relativistic(
        tmp_path, ["particles/relativistic_field_source=mhd"])
    _assert_fatal(result, "relativistic_field_source")


def test_relativistic_hc_rejects_ct_emf_field_source_cpu(tmp_path):
    """The ideal-MHD source does not authorize CT-edge EMF reuse."""
    result = _run_relativistic(
        tmp_path, ["particles/relativistic_field_source=ct_emf"])
    _assert_fatal(result, "relativistic_field_source")


def test_relativistic_hc_coupled_frozen_tn_constructor_runs_without_ce0_cpu(tmp_path):
    """The experimental coupled schema accepts explicit frozen-t^n dynamic ideal MHD."""
    result = _run_coupled(tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr


def test_relativistic_hc_coupled_requires_unit_c_model_cpu(tmp_path):
    """Coupled timesteps remain bounded to the normalized unit-c opening."""
    result = _run_coupled(tmp_path, ["particles/c_model=2.0"])
    _assert_fatal(result, "requires <particles>/c_model = 1.0")


def test_relativistic_hc_coupled_rejects_multiple_same_level_meshblocks_cpu(tmp_path):
    """Same-level migration remains closed until the Phase-8 ownership qualification."""
    result = _run_coupled(tmp_path, ["meshblock/nx1=4"])
    _assert_fatal(result, "requires exactly one meshblock")


def test_relativistic_hc_prescribed_retains_multiple_same_level_meshblocks_cpu(tmp_path):
    """The coupled-only migration fence does not narrow the Phase-4a harness."""
    result = _run_relativistic(
        tmp_path, ["meshblock/nx1=4", "time/nlim=1"])
    assert result.returncode == 0, result.stdout + result.stderr


def test_relativistic_hc_coupled_requires_explicit_temporal_mode_cpu(tmp_path):
    """The coupled grid state is selected explicitly rather than inferred."""
    result = _run_capture(
        tmp_path,
        ["particles/relativistic_field_source=mhd_ideal",
         "time/evolution=dynamic", "mhd/rsolver=hlld"],
        drop_parameters=COUPLED_DROP_PARAMETERS)
    _assert_fatal(result, "relativistic_temporal_sampling = frozen_tn")


def test_relativistic_hc_coupled_rejects_unknown_temporal_mode_cpu(tmp_path):
    """Only the preregistered frozen-t^n coupled temporal model is open."""
    result = _run_coupled(tmp_path, temporal_mode="stage_centered")
    _assert_fatal(result, "relativistic_temporal_sampling")


def test_relativistic_hc_coupled_rejects_obsolete_temporal_spelling_cpu(tmp_path):
    """The preregistered temporal selector has one spelling and no silent alias."""
    result = _run_coupled(
        tmp_path, temporal_mode=None,
        inject_in_particles="relativistic_mhd_temporal_mode = frozen_tn\n")
    _assert_fatal(result, "relativistic_mhd_temporal_mode is unsupported")


def test_relativistic_hc_prescribed_rejects_coupled_temporal_mode_cpu(tmp_path):
    """The coupled temporal selector cannot leak into the prescribed harness."""
    result = _run_relativistic(
        tmp_path, inject_in_particles="relativistic_temporal_sampling = frozen_tn\n")
    _assert_fatal(result, "applies only")


def test_relativistic_hc_coupled_rejects_prescribed_ce0_cpu(tmp_path):
    """Coupled sampled fields are gathered, not seeded from analytical cE constants."""
    result = _run_relativistic(
        tmp_path,
        ["particles/relativistic_field_source=mhd_ideal",
         "time/evolution=dynamic", "mhd/rsolver=hlld"],
        inject_in_particles="relativistic_temporal_sampling = frozen_tn\n")
    _assert_fatal(result, "prescribed-test-only cE0x")


@pytest.mark.parametrize(
    ("extra_flags", "expected_substring"),
    [
        (["time/evolution=kinematic", "mhd/rsolver=advect"], "evolution = dynamic"),
        (["mesh/nx3=1", "meshblock/nx3=1"], "3d"),
        (["mesh/ix1_bc=outflow", "mesh/ox1_bc=outflow"], "periodic"),
        (["particles/interpolation=tsc"], "interpolation"),
    ],
    ids=["dynamic-only", "three-dimensional-only", "periodic-only", "trilinear-only"],
)
def test_relativistic_hc_coupled_rejects_unsupported_modes_cpu(
        tmp_path, extra_flags, expected_substring):
    """Experimental coupled parsing retains the bounded relativistic envelope."""
    result = _run_coupled(tmp_path, extra_flags)
    _assert_fatal(result, expected_substring)


def test_relativistic_hc_coupled_rejects_multilevel_mesh_cpu(tmp_path):
    """The coupled opening remains serial uniform-level only."""
    result = _run_coupled(
        tmp_path,
        append_text=("\n<mesh_refinement>\nrefinement = adaptive\n"
                     "num_levels = 2\nmax_nmb_per_rank = 16\n"))
    _assert_fatal(result, "serial uniform-level mesh")


@pytest.mark.parametrize(
    "append_text",
    [
        "\n<hydro_srcterms>\n",
        "\n<rad_srcterms>\n",
    ],
    ids=["orphan-hydro-source-block", "orphan-radiation-source-block"],
)
def test_relativistic_hc_coupled_rejects_orphan_source_blocks_cpu(
        tmp_path, append_text):
    """Source blocks remain rejected even when their owning physics block is absent."""
    result = _run_coupled(tmp_path, append_text=append_text)
    _assert_fatal(result, "coupled physics blocks")


def test_relativistic_hc_coupled_rejects_user_sources_cpu(tmp_path):
    """Problem-generator user sources cannot bypass the explicit source-block fence."""
    result = _run_coupled(tmp_path, append_text="\n<problem>\nuser_srcs = true\n")
    _assert_fatal(result, "coupled physics blocks")


def test_relativistic_hc_coupled_manufactured_source_exception_is_test_pgen_only_cpu(
        tmp_path):
    """The evolving-field test source is not a general user-source bypass."""
    result = _run_coupled(
        tmp_path,
        append_text=("\n<problem>\nuser_srcs = true\n"
                     "relativistic_coupled_manufactured_test = true\n"))
    _assert_fatal(result, "requires the compiled coupled-runtime test harness")


def test_relativistic_hc_coupled_dependency_is_explicit_cpu():
    """The frozen-t^n task fence is encoded rather than implied by insertion order."""
    source = (ROOT / "src" / "particles" / "particles_tasks.cpp").read_text()
    assert "relativistic_field_source == RelativisticFieldSource::mhd_ideal" in source
    assert "push_dependency = pmy_pack->pmhd->id.savest;" in source
    assert "AddTask(&Particles::Push, this, push_dependency)" in source


def test_relativistic_hc_full_mesh_restart_guard_is_explicit_cpu():
    """Full-mesh -r restoration remains closed until the Phase-6 schema exists."""
    source = (ROOT / "src" / "main.cpp").read_text()
    assert "res_flag && pinput->DoesBlockExist(\"particles\")" in source
    assert "full mesh restart input until the relativistic restart schema is " in source


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


def test_relativistic_hc_positive_cycle_prescribed_test_runs_cpu(tmp_path):
    """The bounded prescribed-test HC pusher is executable after Phase 4a opens."""
    result = _run_relativistic(tmp_path, ["time/nlim=1"])
    assert result.returncode == 0, result.stdout + result.stderr


def test_relativistic_hc_accepts_strict_subcycling_cpu(tmp_path):
    """The Phase-5 acceleration-aware subcycle path is explicit and executable."""
    result = _run_relativistic(tmp_path, inject_in_particles="subcycle = true\n")
    assert result.returncode == 0, result.stdout + result.stderr


def test_relativistic_hc_rejects_subcycle_cap_clipping_cpu(tmp_path):
    """Relativistic subcycling never clips an excessive global request."""
    result = _run_relativistic(
        tmp_path,
        inject_in_particles="subcycle = true\nsubcycle_strict = false\n")
    _assert_fatal(result, "cap clipping is unsupported")


@pytest.mark.parametrize(
    "parameter",
    ["subcycle_momentum_fraction", "subcycle_kinetic_energy_fraction"],
)
def test_relativistic_hc_rejects_deferred_fractional_subcycle_scales_cpu(
        tmp_path, parameter):
    """Near-rest-singular fractional controls remain closed pending a reviewed scale."""
    result = _run_relativistic(
        tmp_path, inject_in_particles=f"{parameter} = 0.1\n")
    _assert_fatal(result, "separately reviewed nonsingular near-rest scale")


def test_relativistic_hc_rejects_multilevel_mesh_cpu(tmp_path):
    """SMR and AMR remain closed until relativistic migration is qualified."""
    result = _run_relativistic(
        tmp_path,
        append_text=("\n<mesh_refinement>\nrefinement = adaptive\n"
                     "num_levels = 2\nmax_nmb_per_rank = 16\n"))
    _assert_fatal(result, "serial uniform-level mesh")


def test_relativistic_hc_rejects_gather_diagnostic_escape_for_normal_pgen_cpu(
        tmp_path):
    """The construction-only Phase-3 escape hatch is not a general runtime opt-out."""
    result = _run_relativistic(
        tmp_path, append_text="\n<problem>\nrelativistic_gather_diagnostic_only = true\n")
    _assert_fatal(result, "requires pgen_name = cr_relativistic_runtime_gather_test")


def test_relativistic_hc_rejects_gather_diagnostic_escape_for_positive_cycle_cpu(
        tmp_path):
    """The registered Phase-3 pgen cannot use the escape hatch during a push cycle."""
    result = _run_relativistic(
        tmp_path,
        ["problem/pgen_name=cr_relativistic_runtime_gather_test", "time/nlim=1"],
        append_text="\n<problem>\nrelativistic_gather_diagnostic_only = true\n")
    _assert_fatal(result, "<time>/nlim = 0")


@pytest.mark.parametrize("b_profile", ["linear_cross", "gradb", "turbulent"])
def test_relativistic_hc_prescribed_test_requires_uniform_b_profile_cpu(
        tmp_path, b_profile):
    """The frozen per-particle prescribed B tuple cannot masquerade as a grid profile."""
    result = _run_relativistic(tmp_path, [f"problem/B_profile={b_profile}"])
    _assert_fatal(result, "B_profile = uniform")


def test_relativistic_hc_rejects_unknown_initial_state_selector_cpu(tmp_path):
    """The authoritative-state selector is an exact closed enum."""
    result = _run_relativistic(
        tmp_path, ["problem/relativistic_initial_state=momentum"])
    _assert_fatal(result, "Unknown relativistic_initial_state")


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("B0x", "nan"),
        ("B0y", "inf"),
        ("B0z", "-inf"),
        ("cE0x", "nan"),
        ("cE0y", "inf"),
        ("cE0z", "-inf"),
    ],
)
def test_relativistic_hc_rejects_nonfinite_prescribed_fields_cpu(
        tmp_path, field, value):
    """Prescribed analytical fields must be finite before device initialization."""
    result = _run_relativistic(tmp_path, [f"problem/{field}={value}"])
    _assert_fatal(result, "fields and initial state must be finite")


@pytest.mark.parametrize("missing_field", ["v0x", "v0y", "v0z"])
def test_relativistic_hc_velocity_initialization_requires_complete_tuple_cpu(
        tmp_path, missing_field):
    """Velocity initialization must be explicit and cannot silently use defaults."""
    result = _run_relativistic(tmp_path, drop_parameters=(missing_field,))
    _assert_fatal(result, "velocity initialization requires")


def test_relativistic_hc_velocity_initialization_rejects_random_mode_cpu(tmp_path):
    """Velocity-mode authoritative state cannot be seeded from random legacy paths."""
    result = _run_relativistic(tmp_path, ["problem/particle_velocity=random"])
    _assert_fatal(result, "particle_velocity = uniform")


def test_relativistic_hc_velocity_initialization_rejects_scalar_v0_cpu(tmp_path):
    """Velocity-mode authoritative state rejects an ambiguous scalar-speed spelling."""
    result = _run_relativistic(tmp_path, append_text="\n<problem>\nv0 = 0.2\n")
    _assert_fatal(result, "no scalar v0")


@pytest.mark.parametrize("partial_w", ["w0x=0.0", "w0y=0.0", "w0z=0.0"])
def test_relativistic_hc_w_initialization_requires_complete_tuple_cpu(
        tmp_path, partial_w):
    """Momentum initialization fails closed on partial authoritative tuples."""
    result = _run_relativistic(
        tmp_path,
        ["problem/relativistic_initial_state=w"],
        drop_parameters=("v0x", "v0y", "v0z"),
        append_text=f"\n<problem>\n{partial_w}\n")
    _assert_fatal(result, "requires exactly w0x, w0y, w0z")


def test_relativistic_hc_w_initialization_runs_cpu(tmp_path):
    """Explicit authoritative momentum initialization reaches the bounded HC path."""
    result = _run_relativistic(
        tmp_path,
        ["problem/relativistic_initial_state=w", "time/nlim=1"],
        drop_parameters=("v0x", "v0y", "v0z"),
        append_text="\n<problem>\nw0x = 0.2\nw0y = -0.1\nw0z = 0.05\n")
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("value", ["nan", "inf", "-inf"])
def test_relativistic_hc_rejects_nonfinite_w_initialization_cpu(tmp_path, value):
    """Authoritative momentum tuples must be finite and shadow-representable."""
    result = _run_relativistic(
        tmp_path,
        ["problem/relativistic_initial_state=w"],
        drop_parameters=("v0x", "v0y", "v0z"),
        append_text=f"\n<problem>\nw0x = {value}\nw0y = 0.0\nw0z = 0.0\n")
    _assert_fatal(result, "fields and initial state must be finite")


def test_relativistic_hc_w_initialization_rejects_velocity_tuple_cpu(tmp_path):
    """Momentum mode cannot retain a contradictory legacy velocity tuple."""
    result = _run_relativistic(
        tmp_path,
        ["problem/relativistic_initial_state=w"],
        append_text="\n<problem>\nw0x = 0.2\nw0y = 0.0\nw0z = 0.0\n")
    _assert_fatal(result, "and no velocity initialization values")


def test_relativistic_hc_velocity_initialization_rejects_w_tuple_cpu(tmp_path):
    """Velocity mode cannot retain a contradictory authoritative momentum tuple."""
    result = _run_relativistic(
        tmp_path, append_text="\n<problem>\nw0x = 0.2\nw0y = 0.0\nw0z = 0.0\n")
    _assert_fatal(result, "velocity initialization rejects")


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


def test_legacy_pushers_reject_relativistic_temporal_sampling_cpu(tmp_path):
    """Legacy modes must reject the solver-coupled temporal selector."""
    flags = [
        "particles/pusher=boris",
        "time/evolution=dynamic",
        "mhd/rsolver=hlld",
        "time/nlim=1",
    ]
    result = _run_capture(
        tmp_path, flags,
        input_file=LEGACY_INPUT,
        inject_in_particles="relativistic_temporal_sampling = frozen_tn\n")
    _assert_fatal(result, "relativistic-only")


@pytest.mark.parametrize(
    "problem_key",
    [
        "relativistic_initial_state = velocity\n",
        "w0x = 0.0\n",
        "w0y = 0.0\n",
        "w0z = 0.0\n",
        "cE0x = 0.0\n",
        "cE0y = 0.0\n",
        "cE0z = 0.0\n",
        "relativistic_gather_diagnostic_only = true\n",
    ],
)
def test_legacy_pushers_reject_relativistic_problem_keys_cpu(tmp_path, problem_key):
    """Legacy modes must reject Phase-4a prescribed-test initialization keys."""
    flags = [
        "particles/pusher=boris",
        "time/evolution=dynamic",
        "mhd/rsolver=hlld",
        "time/nlim=1",
    ]
    result = _run_capture(tmp_path, flags, append_text=f"\n<problem>\n{problem_key}")
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
