"""
Regression tests for multigrid self-gravity on serial CPU builds.

These tests cover the supported user-facing contracts: periodic Jeans solves,
isolated multipole solves, source-term coupling, restart continuity, AMR smoke,
profile diagnostics, and explicit failures for invalid inputs.
"""

from pathlib import Path
import math

import numpy as np
import pytest

import test_suite.testutils as testutils
from test_suite.selfgravity import selfgravity_utils as sg


def _remove_block(text, block_name):
    lines = text.splitlines()
    output = []
    skipping = False
    for line in lines:
        stripped = line.strip()
        if stripped == f"<{block_name}>":
            skipping = True
            continue
        if skipping and stripped.startswith("<") and stripped.endswith(">"):
            skipping = False
        if not skipping:
            output.append(line)
    return "\n".join(output) + "\n"


def _remove_parameters(text, block_name, parameters):
    lines = text.splitlines()
    output = []
    in_block = False
    for line in lines:
        stripped = line.strip()
        if stripped == f"<{block_name}>":
            in_block = True
            output.append(line)
            continue
        if in_block and stripped.startswith("<") and stripped.endswith(">"):
            in_block = False
        key = stripped.split("=", 1)[0].strip()
        if in_block and key in parameters:
            continue
        output.append(line)
    return "\n".join(output) + "\n"


def _write_transformed_input(tmp_path, source, transform):
    text = (sg.REPO_ROOT / source).read_text()
    path = tmp_path / Path(source).name
    path.write_text(transform(text))
    return path


@pytest.mark.parametrize(
    "input_file,basename",
    [
        ("inputs/tests/selfgravity.athinput", "selfgravity_hydro_cpu"),
        ("inputs/tests/selfgravity_mhd.athinput", "selfgravity_mhd_cpu"),
    ],
)
def test_jeans_wave_selfgravity_converges(input_file, basename):
    try:
        sg.run_athena(input_file, basename, expect_jeans_banner=True)
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_periodic_jeans_potential_matches_analytic_solution():
    basename = "selfgravity_phi_analytic_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity.athinput", basename, expect_jeans_banner=True
        )
        output = sg.parse_binary_output(sg.latest_binary_output(basename, "grav_phi"))
        phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")

        k_wave, four_pi_g, _ = sg.jeans_geometry(n_jeans=0.5)
        phase = sg.jeans_phase(x1, x2, x3)
        analytic = -four_pi_g * 1.0e-3 / (k_wave * k_wave) * np.sin(k_wave * phase)

        numerical = phi - np.mean(phi)
        analytic -= np.mean(analytic)
        assert sg.relative_l2(numerical - analytic, analytic) < 3.0e-2
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_initial_output_contains_solved_jeans_potential():
    basename = "selfgravity_phi_initial_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity.athinput", basename, expect_jeans_banner=True
        )
        output = sg.parse_binary_output(sg.first_binary_output(basename, "grav_phi"))
        phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")

        k_wave, four_pi_g, _ = sg.jeans_geometry(n_jeans=0.5)
        phase = sg.jeans_phase(x1, x2, x3)
        analytic = -four_pi_g * 1.0e-3 / (k_wave * k_wave) * np.sin(k_wave * phase)

        numerical = phi - np.mean(phi)
        analytic -= np.mean(analytic)
        assert output["metadata"]["cycle"] == 0
        assert sg.relative_l2(numerical - analytic, analytic) < 3.0e-2
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_final_output_potential_matches_final_density():
    basename = "selfgravity_phi_final_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            basename,
            extra_args=["time/cfl_number=0.2", "problem/amp=1e-3"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        phi_out = sg.parse_binary_output(sg.latest_binary_output(basename, "grav_phi"))
        hydro_out = sg.parse_binary_output(sg.latest_binary_output(basename, "hydro_u"))
        phi, x1, x2, x3 = sg.uniform_grid_field(phi_out, "grav_phi")
        dens, _, _, _ = sg.uniform_grid_field(hydro_out, "dens")
        dx1 = x1[1] - x1[0]
        dx2 = x2[1] - x2[0]
        dx3 = x3[1] - x3[0]
        laplacian = (
            (np.roll(phi, -1, axis=2) - 2.0 * phi + np.roll(phi, 1, axis=2)) / (dx1 * dx1)
            + (np.roll(phi, -1, axis=1) - 2.0 * phi + np.roll(phi, 1, axis=1))
            / (dx2 * dx2)
            + (np.roll(phi, -1, axis=0) - 2.0 * phi + np.roll(phi, 1, axis=0))
            / (dx3 * dx3)
        )
        _, four_pi_g, _ = sg.jeans_geometry(n_jeans=0.5)
        source = four_pi_g * (dens - np.mean(dens))
        assert phi_out["metadata"]["cycle"] == 1
        assert hydro_out["metadata"]["cycle"] == 1
        assert sg.relative_l2(laplacian - source, source) < 1.0e-4
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_zero_field_mhd_matches_hydro_gravity_solution():
    hydro_basename = "selfgravity_equiv_hydro_cpu"
    mhd_basename = "selfgravity_equiv_mhd_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity.athinput", hydro_basename, expect_jeans_banner=True
        )
        sg.run_athena(
            "inputs/tests/selfgravity_mhd.athinput",
            mhd_basename,
            expect_jeans_banner=True,
        )
        hydro_output = sg.parse_binary_output(
            sg.latest_binary_output(hydro_basename, "grav_phi")
        )
        mhd_output = sg.parse_binary_output(
            sg.latest_binary_output(mhd_basename, "grav_phi")
        )
        hydro_phi, _, _, _ = sg.concatenate_field(hydro_output, "grav_phi")
        mhd_phi, _, _, _ = sg.concatenate_field(mhd_output, "grav_phi")
        assert np.max(np.abs(hydro_phi - mhd_phi)) < 1.0e-10
    finally:
        sg.cleanup_outputs(hydro_basename, mhd_basename)
        testutils.cleanup()


def test_source_terms_follow_analytic_jeans_force():
    on_basename = "selfgravity_source_on_cpu"
    off_basename = "selfgravity_source_off_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            on_basename,
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        sg.run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            off_basename,
            extra_args=["hydro_srcterms/self_gravity=false"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        on = sg.parse_binary_output(sg.latest_binary_output(on_basename, "hydro_u"))
        off = sg.parse_binary_output(sg.latest_binary_output(off_basename, "hydro_u"))
        time = on["metadata"]["time"]
        dens, x1, x2, x3 = sg.concatenate_field(on, "dens")
        k_wave, four_pi_g, direction = sg.jeans_geometry(n_jeans=0.5)
        phase = sg.jeans_phase(x1, x2, x3)
        accel = (four_pi_g * 1.0e-4 / k_wave) * np.cos(k_wave * phase)
        expected = [dens * accel * direction[n] * time for n in range(3)]

        for name, ref in zip(("mom1", "mom2", "mom3"), expected):
            mom_on, _, _, _ = sg.concatenate_field(on, name)
            mom_off, _, _, _ = sg.concatenate_field(off, name)
            assert sg.relative_l2((mom_on - mom_off) - ref, ref) < 1.5e-1
    finally:
        sg.cleanup_outputs(on_basename, off_basename)
        testutils.cleanup()


def _uniform_sphere_potential(x, y, z, center, mass, radius, gconst):
    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]
    r = np.sqrt(dx * dx + dy * dy + dz * dz)
    outside = -gconst * mass / np.maximum(r, 1.0e-300)
    inside = -gconst * mass * (3.0 * radius * radius - r * r) / (2.0 * radius**3)
    return np.where(r < radius, inside, outside)


def _point_mass_acceleration(x, y, z, center, mass, gconst):
    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]
    r2 = dx * dx + dy * dy + dz * dz
    r3 = np.maximum(r2, 1.0e-300) ** 1.5
    return (-gconst * mass * dx / r3, -gconst * mass * dy / r3, -gconst * mass * dz / r3)


def test_binary_multipole_potential_and_force_are_physical():
    basename = "selfgravity_binary_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity_binary.athinput", basename, max_final_defect=2.0e-6
        )
        output = sg.parse_binary_output(sg.latest_binary_output(basename, "grav_phi"))
        phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")
        gconst = 1.0 / (4.0 * math.pi)
        radius = 0.10
        center1 = (0.15, 0.0, 0.0)
        center2 = (-0.15, 0.0, 0.0)
        analytic = _uniform_sphere_potential(
            x1, x2, x3, center1, 1.0, radius, gconst
        ) + _uniform_sphere_potential(x1, x2, x3, center2, 0.5, radius, gconst)
        r1 = np.sqrt((x1 - center1[0]) ** 2 + x2 * x2 + x3 * x3)
        r2 = np.sqrt((x1 - center2[0]) ** 2 + x2 * x2 + x3 * x3)
        far_source = (r1 > 1.5 * radius) & (r2 > 1.5 * radius)
        assert (
            sg.relative_l2(phi[far_source] - analytic[far_source], analytic[far_source])
            < 3.5e-1
        )

        phi_grid, gx1, gx2, gx3 = sg.uniform_grid_field(output, "grav_phi")
        dphidz, dphidy, dphidx = np.gradient(phi_grid, gx3, gx2, gx1, edge_order=2)
        x3g, x2g, x1g = np.meshgrid(gx3, gx2, gx1, indexing="ij")
        ax1a, ax2a, ax3a = _point_mass_acceleration(x1g, x2g, x3g, center1, 1.0, gconst)
        bx1a, bx2a, bx3a = _point_mass_acceleration(x1g, x2g, x3g, center2, 0.5, gconst)
        analytic_ax = ax1a + bx1a
        analytic_ay = ax2a + bx2a
        analytic_az = ax3a + bx3a
        r1_grid = np.sqrt((x1g - center1[0]) ** 2 + x2g * x2g + x3g * x3g)
        r2_grid = np.sqrt((x1g - center2[0]) ** 2 + x2g * x2g + x3g * x3g)
        mask = (r1_grid > 2.0 * radius) & (r2_grid > 2.0 * radius)
        mask[:, :, :1] = False
        mask[:, :, -1:] = False
        mask[:, :1, :] = False
        mask[:, -1:, :] = False
        mask[:1, :, :] = False
        mask[-1:, :, :] = False
        numerical = np.stack((-dphidx[mask], -dphidy[mask], -dphidz[mask]))
        reference = np.stack((analytic_ax[mask], analytic_ay[mask], analytic_az[mask]))
        cosine = np.sum(numerical * reference) / (
            np.linalg.norm(numerical) * np.linalg.norm(reference)
        )
        assert cosine > 0.90
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_restart_equivalence_for_self_gravity_cycle():
    straight = "selfgravity_restart_straight_cpu"
    split = "selfgravity_restart_split_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            straight,
            extra_args=["time/nlim=2"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        sg.run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            split,
            extra_args=["time/nlim=1"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        restart = sg.latest_restart_output(split)
        sg.run_restart(restart, extra_args=["time/nlim=2"])

        straight_out = sg.parse_binary_output(
            sg.latest_binary_output(straight, "hydro_u")
        )
        split_out = sg.parse_binary_output(sg.latest_binary_output(split, "hydro_u"))
        for label in ("dens", "mom1", "mom2", "mom3"):
            straight_field, _, _, _ = sg.concatenate_field(straight_out, label)
            split_field, _, _, _ = sg.concatenate_field(split_out, label)
            assert np.max(np.abs(straight_field - split_field)) < 1.0e-6
    finally:
        sg.cleanup_outputs(straight, split)
        testutils.cleanup()


def test_root_grid_can_run_on_host():
    basename = "selfgravity_root_host_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity.athinput",
            basename,
            extra_args=["gravity/root_on_host=true"],
            expect_jeans_banner=True,
        )
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_fixed_iteration_mode_runs():
    basename = "selfgravity_fixed_iter_cpu"
    try:
        sg.run_athena(
            "inputs/tests/selfgravity.athinput",
            basename,
            extra_args=["gravity/threshold=-1", "gravity/niteration=4"],
            max_final_defect=None,
            expect_jeans_banner=True,
        )
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_profile_diagnostics_report_gravity_phases():
    basename = "selfgravity_profile_cpu"
    try:
        output = sg.run_athena(
            "inputs/tests/selfgravity.athinput",
            basename,
            extra_args=["gravity/profile=true"],
            expect_jeans_banner=True,
        )
        assert "[gravity-profile]" in output
        profile_lines = [
            line for line in output.splitlines() if line.startswith("[gravity-profile]")
        ]
        assert [line.split("stage=", 1)[1].split()[0] for line in profile_lines] == [
            "0",
            "2",
            "3",
        ]
        for field in (
            "source_load=",
            "setup=",
            "root_transfer=",
            "smooth=",
            "boundary=",
            "restrict_prolong=",
            "result_retrieve=",
        ):
            assert field in output
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_be_collapse_amr_smoke():
    basename = "selfgravity_be_amr_cpu"
    try:
        output = sg.run_athena(
            "inputs/tests/selfgravity_be_amr.athinput",
            basename,
            max_final_defect=2.0e-5,
        )
        assert "Number of physical levels of refinement = 1" in output
        parsed = sg.parse_binary_output(sg.latest_binary_output(basename, "grav_phi"))
        assert parsed["blocks"]
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


def test_dynamic_amr_regrid_smoke(tmp_path):
    def append_adaptive_blocks(text):
        return (
            text.rstrip()
            + """

<mesh_refinement>
refinement = adaptive
max_nmb_per_rank = 128
num_levels = 2
ncycle_check = 1
refinement_interval = 1

<amr_criterion0>
method = min_max
variable = hydro_w_d
value_max = 0.5
"""
        )

    basename = "selfgravity_dynamic_amr_cpu"
    path = _write_transformed_input(
        tmp_path,
        "inputs/tests/selfgravity.athinput",
        append_adaptive_blocks,
    )
    try:
        output = sg.run_athena(
            path,
            basename,
            extra_args=["mesh/nghost=4", "time/nlim=2", "gravity/threshold=1e-6"],
            max_final_defect=2.0e-6,
            expect_jeans_banner=True,
        )
        assert "Current number of MeshBlocks = 64" in output
        assert "56 MeshBlocks created, 0 deleted by AMR" in output
        parsed = sg.parse_binary_output(sg.latest_binary_output(basename, "grav_phi"))
        assert len(parsed["blocks"]) == 64
    finally:
        sg.cleanup_outputs(basename)
        testutils.cleanup()


@pytest.mark.parametrize(
    "input_file,extra_args,expected",
    [
        (
            "inputs/tests/selfgravity.athinput",
            ["meshblock/nx3=8"],
            "logically cubic MeshBlocks",
        ),
        (
            "inputs/tests/selfgravity.athinput",
            ["mesh/nx3=1", "meshblock/nx3=1"],
            "requires a 3D mesh",
        ),
        (
            "inputs/tests/selfgravity.athinput",
            ["gravity/mg_bc=bogus"],
            "Unknown gravity/mg_bc",
        ),
        (
            "inputs/tests/selfgravity_binary.athinput",
            ["gravity/mporder=3"],
            "gravity/mporder must be 2 or 4",
        ),
    ],
)
def test_invalid_runtime_input_fails(input_file, extra_args, expected):
    try:
        sg.expect_failure(input_file, extra_args=extra_args, expected_text=expected)
    finally:
        testutils.cleanup()


def test_missing_gravity_block_fails(tmp_path):
    path = _write_transformed_input(
        tmp_path,
        "inputs/tests/selfgravity.athinput",
        lambda text: _remove_block(text, "gravity"),
    )
    try:
        sg.expect_failure(path, expected_text="self_gravity source terms require")
    finally:
        testutils.cleanup()


def test_missing_convergence_control_fails(tmp_path):
    path = _write_transformed_input(
        tmp_path,
        "inputs/tests/selfgravity.athinput",
        lambda text: _remove_parameters(text, "gravity", {"threshold", "niteration"}),
    )
    try:
        sg.expect_failure(
            path, expected_text="Set either gravity/threshold or gravity/niteration"
        )
    finally:
        testutils.cleanup()


def test_two_fluid_self_gravity_fails(tmp_path):
    def append_two_fluid_blocks(text):
        return (
            text.rstrip()
            + """

<hydro>
eos = isothermal
reconstruct = plm
rsolver = llf
iso_sound_speed = 1.0

<ion-neutral>
drag_coeff = 1.0
"""
        )

    path = _write_transformed_input(
        tmp_path,
        "inputs/tests/selfgravity_mhd.athinput",
        append_two_fluid_blocks,
    )
    try:
        sg.expect_failure(path, expected_text="Self-gravity currently supports one fluid")
    finally:
        testutils.cleanup()
