#!/usr/bin/env python3
"""Run and summarize reproducible CGL Landau-fluid validation workflows."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT_DIR = Path(__file__).resolve().parents[1]
SCHEMA_VERSION = 2
HISTORY_LABEL = re.compile(r"\[\d+\]=(\S+)")
FULL_INPUTS = (
    "inputs/unit_tests/cgl_lf_quant_parallel.athinput",
    "inputs/unit_tests/cgl_lf_quant_parallel_collisional.athinput",
    "inputs/unit_tests/cgl_lf_quant_perp.athinput",
    "inputs/unit_tests/cgl_lf_quant_perp_collisional.athinput",
    "inputs/unit_tests/cgl_lf_quant_grad_b.athinput",
    "inputs/unit_tests/cgl_lf_flux_limiter.athinput",
    "inputs/unit_tests/cgl_lf_limiter_heat_flux_suppression.athinput",
    "inputs/unit_tests/cgl_lf_limiter_mirror.athinput",
    "inputs/unit_tests/cgl_lf_limiter_firehose.athinput",
    "inputs/unit_tests/cgl_lf_field_wave.athinput",
    "inputs/unit_tests/cgl_lf_paper_oblique_wave.athinput",
    "inputs/unit_tests/cgl_pure_paper_oblique_wave.athinput",
    "inputs/unit_tests/cgl_pure_paper_eigen_alfven.athinput",
    "inputs/unit_tests/cgl_pure_paper_eigen_slow.athinput",
    "inputs/unit_tests/cgl_pure_paper_eigen_fast.athinput",
    "inputs/unit_tests/cgl_lf_paper_eigen_alfven.athinput",
    "inputs/unit_tests/cgl_lf_paper_eigen_slow.athinput",
    "inputs/unit_tests/cgl_lf_paper_eigen_fast.athinput",
)
UNSAFE_LF_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
)


@dataclass(frozen=True)
class CaseSpec:
    """Description of one AthenaK run within a workflow."""

    name: str
    input_path: str
    overrides: tuple[str, ...] = ()
    lf_active: bool = True
    validation_output: bool = False
    comparison_group: str | None = None
    comparison_role: str | None = None
    amr: bool = False
    paper_smoke: bool = False
    accuracy_study: str | None = None
    accuracy_parameters: tuple[tuple[str, str], ...] = ()


@dataclass
class RunPaths:
    """Paths owned by one workflow result bundle."""

    root: Path
    logs: Path = field(init=False)
    history: Path = field(init=False)
    data: Path = field(init=False)
    figures: Path = field(init=False)
    inputs: Path = field(init=False)
    cases: Path = field(init=False)

    def __post_init__(self) -> None:
        self.logs = self.root / "logs"
        self.history = self.root / "history"
        self.data = self.root / "data"
        self.figures = self.root / "figures"
        self.inputs = self.root / "inputs"
        self.cases = self.root / "cases"

    def create(self) -> None:
        for path in (
            self.root,
            self.logs,
            self.history,
            self.data,
            self.figures,
            self.inputs,
            self.cases,
        ):
            path.mkdir(parents=True, exist_ok=True)


def display_path(path: Path, base: Path) -> str:
    """Return a readable relative path when path is inside base."""

    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def resolve_from_root(value: str | Path) -> Path:
    """Resolve relative command-line paths from the repository root."""

    path = Path(value)
    return path if path.is_absolute() else ROOT_DIR / path


def timestamp() -> str:
    """UTC timestamp appropriate for run-directory names and manifests."""

    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def parse_history(path: Path) -> dict[str, list[float]]:
    """Parse an AthenaK history file without importing analysis packages."""

    labels: list[str] = []
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#  [1]="):
                labels = HISTORY_LABEL.findall(line)
            elif line and not line.startswith("#") and line.strip():
                rows.append([float(value) for value in line.split()])
    if not labels or not rows:
        raise ValueError(f"history output is incomplete: {path}")
    return {
        label: [row[index] for row in rows]
        for index, label in enumerate(labels)
    }


def parse_table(path: Path) -> dict[str, list[float]]:
    """Parse an AthenaK formatted-table output."""

    labels: list[str] = []
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# gid"):
                labels = line[1:].split()
            elif line and not line.startswith("#") and line.strip():
                rows.append([float(value) for value in line.split()])
    if not labels or not rows:
        raise ValueError(f"table output is incomplete: {path}")
    return {
        label: [row[index] for row in rows]
        for index, label in enumerate(labels)
    }


def final_value(history: dict[str, list[float]], key: str) -> float:
    """Return the last value of a named history column."""

    if key not in history:
        raise ValueError(f"required history column is missing: {key}")
    return history[key][-1]


def lf_diagnostics(path: Path) -> dict[str, object]:
    """Extract LF safety and limiter occupancy diagnostics from history."""

    history = parse_history(path)
    nstage = final_value(history, "lf_nstage")
    counters = {
        key: final_value(history, key)
        for key in (
            "lf_nstage",
            "lf_dfloor",
            "lf_pfloor",
            "lf_nonfin",
            "lf_nonpos",
            "lf_mirror",
            "lf_firehs",
            "lf_hardbd",
        )
    }
    for key in (
        "lf_qface", "lf_qprcap", "lf_qpr10", "lf_qpecap", "lf_qpe10",
        "lf_qprwrk", "lf_qpewrk", "lf_cpwrk", "lf_cawrk",
    ):
        if key in history:
            counters[key] = final_value(history, key)
    clean = all(counters[key] == 0.0 for key in UNSAFE_LF_COLUMNS)
    denominator = max(nstage, 1.0)
    result: dict[str, object] = {
        "clean": clean,
        "counters": counters,
        "mirror_occupancy": counters["lf_mirror"] / denominator,
        "firehose_occupancy": counters["lf_firehs"] / denominator,
    }
    if counters.get("lf_qface", 0.0) > 0.0:
        qfaces = counters["lf_qface"]
        result["heat_flux_cap_fractions"] = {
            "parallel_over_1": counters["lf_qprcap"] / qfaces,
            "parallel_over_10": counters["lf_qpr10"] / qfaces,
            "perpendicular_over_1": counters["lf_qpecap"] / qfaces,
            "perpendicular_over_10": counters["lf_qpe10"] / qfaces,
        }
    if "lf_qprwrk" in counters and "lf_qpewrk" in counters:
        result["heat_flux_work"] = {
            "definition": (
                "Cumulative signed RKL2-applied owned-face contraction "
                "of capped parallel/perpendicular LF heat fluxes"
            ),
            "amr_face_ownership": "coarse/fine interfaces use the fine-side flux",
            "signed": True,
            "parallel": counters["lf_qprwrk"],
            "perpendicular": counters["lf_qpewrk"],
            "total": counters["lf_qprwrk"] + counters["lf_qpewrk"],
        }
    if "lf_cpwrk" in counters and "lf_cawrk" in counters:
        result["pressure_work"] = {
            "definition": (
                "Cumulative RK-applied cell contraction of velocity with the "
                "AMR-corrected CGL pressure-traction momentum flux divergence"
            ),
            "anisotropic_definition": "retained Delta-p traction contribution",
            "signed": True,
            "total": counters["lf_cpwrk"],
            "anisotropic": counters["lf_cawrk"],
        }
    return result


def append_history_output(source_text: str) -> str:
    """Ensure a staged input emits history diagnostics."""

    has_hst = re.search(r"(?ms)^<output\d+>.*?^\s*file_type\s*=\s*hst\b",
                        source_text)
    if has_hst:
        return source_text
    output_ids = [
        int(value)
        for value in re.findall(r"(?m)^<output(\d+)>", source_text)
    ]
    output_id = max(output_ids, default=0) + 1
    return (
        source_text.rstrip()
        + "\n\n"
        + f"<output{output_id}>\n"
        + "file_type = hst\n"
        + "data_format = %24.16e\n"
        + "dcycle = 1\n"
    )


def input_block_value(source_text: str, block: str, parameter: str) -> str | None:
    """Read one simple Athena input value from a named block."""

    match = re.search(
        rf"(?ms)^<{re.escape(block)}>\s*(.*?)(?=^<[^>]+>|\Z)",
        source_text,
    )
    if match is None:
        return None
    value = re.search(
        rf"(?m)^\s*{re.escape(parameter)}\s*=\s*([^\s#]+)",
        match.group(1),
    )
    return value.group(1) if value is not None else None


def model_choices(source_text: str, overrides: list[str]) -> dict[str, str]:
    """Archive policy choices that affect interpretation of LF diagnostics."""

    def choice(block: str, parameter: str, default: str) -> str:
        value = input_block_value(source_text, block, parameter) or default
        prefix = f"{block}/{parameter}="
        for override in overrides:
            if override.startswith(prefix):
                value = override[len(prefix):]
        return value

    choices = {
        "passive_delta": choice("mhd", "passive", "false"),
        "paper_mode": choice("problem", "paper_mode", "unspecified"),
        "beta0": choice("problem", "beta0", "unspecified"),
        "mesh_nx1": choice("mesh", "nx1", "unspecified"),
        "mesh_nx2": choice("mesh", "nx2", "unspecified"),
        "mesh_nx3": choice("mesh", "nx3", "unspecified"),
        "domain_x1": (
            f"{choice('mesh', 'x1min', 'unspecified')}:"
            f"{choice('mesh', 'x1max', 'unspecified')}"
        ),
        "domain_x2": (
            f"{choice('mesh', 'x2min', 'unspecified')}:"
            f"{choice('mesh', 'x2max', 'unspecified')}"
        ),
        "domain_x3": (
            f"{choice('mesh', 'x3min', 'unspecified')}:"
            f"{choice('mesh', 'x3max', 'unspecified')}"
        ),
        "analysis_t_start": choice("problem", "analysis_t_start", "unspecified"),
        "analysis_t_end": choice("problem", "analysis_t_end", "unspecified"),
        "time_tlim": choice("time", "tlim", "unspecified"),
        "output1_file_type": choice("output1", "file_type", "unspecified"),
        "output1_dt": choice("output1", "dt", "unspecified"),
        "output1_dcycle": choice("output1", "dcycle", "unspecified"),
        "output2_file_type": choice("output2", "file_type", "unspecified"),
        "output2_dt": choice("output2", "dt", "unspecified"),
        "output2_dcycle": choice("output2", "dcycle", "unspecified"),
        "output2_single_file_per_rank": choice(
            "output2", "single_file_per_rank", "false"
        ),
        "output3_file_type": choice("output3", "file_type", "unspecified"),
        "output3_dt": choice("output3", "dt", "unspecified"),
        "output3_dcycle": choice("output3", "dcycle", "unspecified"),
        "output3_single_file_per_rank": choice(
            "output3", "single_file_per_rank", "false"
        ),
        "cgl_firehose_threshold": choice("mhd", "cgl_firehose_threshold", "oblique"),
        "cgl_heat_flux_integrator": choice("mhd", "cgl_heat_flux_integrator", "sts"),
        "cgl_lf_record_pressure_work": choice(
            "mhd", "cgl_lf_record_pressure_work", "false"
        ),
        "lf_k_parallel": choice("mhd", "lf_k_parallel", "unspecified"),
        "lf_coefficient_mode": choice("mhd", "lf_coefficient_mode", "local"),
        "lf_c_parallel0": choice("mhd", "lf_c_parallel0", "unspecified"),
        "nu_coll": choice("mhd", "nu_coll", "0.0"),
        "mirror_limiter": choice("mhd", "mirror_limiter", "false"),
        "firehose_limiter": choice("mhd", "firehose_limiter", "false"),
        "limiter_nu_coll": choice("mhd", "limiter_nu_coll", "0.0"),
        "backup_limiters": choice("mhd", "backup_limiters", "false"),
        "dfloor": choice("mhd", "dfloor", "unspecified"),
        "pfloor": choice("mhd", "pfloor", "unspecified"),
        "tfloor": choice("mhd", "tfloor", "unspecified"),
        "bfloor": choice("mhd", "bfloor", "unspecified"),
        "cgl_collision_split": "two_half_steps_after_lf_sweeps",
    }
    driving_type = input_block_value(source_text, "turb_driving", "driving_type")
    if driving_type is not None:
        driving_type = choice("turb_driving", "driving_type", driving_type)
        choices.update({
            "forcing_mode": (
                "alfvenic_z_perpendicular" if driving_type == "1"
                else "isotropic_random" if driving_type == "0"
                else f"unsupported_{driving_type}"
            ),
            "forcing_seed": choice("turb_driving", "rseed", "-1"),
            "forcing_dedt": choice("turb_driving", "dedt", "0.0"),
            "forcing_tcorr": choice("turb_driving", "tcorr", "0.0"),
            "forcing_record_injected_work": choice(
                "turb_driving", "record_injected_work", "false"
            ),
            "forcing_nlow": choice("turb_driving", "nlow", "1"),
            "forcing_nhigh": choice("turb_driving", "nhigh", "2"),
            "forcing_physical_k_shell": choice(
                "turb_driving", "physical_k_shell", "false"
            ),
            "forcing_k_shell_unit": choice("turb_driving", "k_shell_unit", "0.0"),
            "forcing_isotropic_power_spectrum": choice(
                "turb_driving", "isotropic_power_spectrum", "false"
            ),
            "forcing_expo": choice("turb_driving", "expo", "5.0/3.0"),
            "forcing_exp_prp": choice("turb_driving", "exp_prp", "5.0/3.0"),
            "forcing_exp_prl": choice("turb_driving", "exp_prl", "0.0"),
            "forcing_time_integration": "ou_once_per_cycle_rk_companion_source_work",
        })
    return choices


def stage_input(spec: CaseSpec, paths: RunPaths) -> Path:
    """Create the execution input retained with a result bundle."""

    source = ROOT_DIR / spec.input_path
    staged = paths.inputs / f"{spec.name}.athinput"
    staged.write_text(
        append_history_output(source.read_text(encoding="utf-8")),
        encoding="utf-8",
    )
    return staged


def workflow_cases(workflow: str) -> list[CaseSpec]:
    """Return the concrete cases for an executable workflow."""

    if workflow == "quick":
        return [
            CaseSpec(
                "cgl_lf_quick",
                "inputs/unit_tests/cgl_lf_quant_parallel.athinput",
                validation_output=True,
            )
        ]
    if workflow == "amr":
        return [
            CaseSpec(
                "cgl_lf_amr",
                "inputs/tests/cgl_lf_amr_2d.athinput",
                amr=True,
            )
        ]
    if workflow == "paper-smoke":
        source = "inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput"
        passive_source = "inputs/cgl_lf_paper/cgl_lf_paper_smoke_passive_beta10.athinput"
        return [
            CaseSpec(
                "paper_smoke_active_alfvenic",
                source,
                paper_smoke=True,
            ),
            CaseSpec(
                "paper_smoke_active_random",
                source,
                (
                    "turb_driving/driving_type=0",
                    "turb_driving/rseed=314159",
                    "turb_driving/expo=2.0",
                ),
                paper_smoke=True,
            ),
            CaseSpec(
                "paper_smoke_passive_alfvenic",
                passive_source,
                paper_smoke=True,
            ),
        ]
    if workflow == "paper-standard":
        prefix = "inputs/cgl_lf_paper/"
        return [
            CaseSpec(
                "paper_standard_active_alfvenic_beta1",
                prefix + "cgl_lf_paper_standard_active_alfvenic_beta1.athinput",
            ),
            CaseSpec(
                "paper_standard_active_alfvenic_beta10",
                prefix + "cgl_lf_paper_standard_active_alfvenic_beta10.athinput",
            ),
            CaseSpec(
                "paper_standard_active_alfvenic_beta100",
                prefix + "cgl_lf_paper_standard_active_alfvenic_beta100.athinput",
            ),
            CaseSpec(
                "paper_standard_passive_alfvenic_beta10",
                prefix + "cgl_lf_paper_standard_passive_alfvenic_beta10.athinput",
            ),
            CaseSpec(
                "paper_standard_active_random_beta10",
                prefix + "cgl_lf_paper_standard_active_random_beta10.athinput",
            ),
        ]
    if workflow == "paper-convergence":
        source = (
            "inputs/cgl_lf_paper/"
            "cgl_lf_paper_standard_active_alfvenic_beta10.athinput"
        )

        def reduced_overrides(resolution: int, *extra: str) -> tuple[str, ...]:
            return (
                f"mesh/nx1={resolution}",
                f"mesh/nx2={resolution}",
                f"mesh/nx3={2*resolution}",
                f"meshblock/nx1={resolution}",
                f"meshblock/nx2={resolution}",
                f"meshblock/nx3={2*resolution}",
                "time/tlim=0.02",
                "time/ndiag=1",
                "problem/analysis_t_start=0.0",
                "problem/analysis_t_end=0.02",
                "output1/dt=0.005",
                "output2/dt=0.01",
                "output3/dt=0.02",
                *extra,
            )

        return [
            CaseSpec(
                f"paper_convergence_resolution_{resolution}",
                source,
                reduced_overrides(resolution),
                paper_smoke=True,
            )
            for resolution in (8, 16, 32)
        ] + [
            CaseSpec(
                "paper_convergence_heat_flux_strong",
                source,
                reduced_overrides(16, "mhd/lf_k_parallel=0.06283185307179586"),
                paper_smoke=True,
            ),
            CaseSpec(
                "paper_convergence_heat_flux_weak",
                source,
                reduced_overrides(16, "mhd/lf_k_parallel=628.3185307179587"),
                paper_smoke=True,
            ),
            CaseSpec(
                "paper_convergence_firehose_oblique",
                source,
                reduced_overrides(16, "mhd/cgl_firehose_threshold=oblique"),
                paper_smoke=True,
            ),
        ]
    if workflow == "paper-nulim":
        prefix = "inputs/cgl_lf_paper/"
        return [
            CaseSpec(
                "paper_nulim_beta100_20",
                prefix + "cgl_lf_paper_nulim_beta100_20.athinput",
            ),
            CaseSpec(
                "paper_nulim_beta100_200",
                prefix + "cgl_lf_paper_nulim_beta100_200.athinput",
            ),
            CaseSpec(
                "paper_nulim_beta100_hardwall",
                prefix + "cgl_lf_paper_nulim_beta100_hardwall.athinput",
            ),
        ]
    if workflow == "accuracy":
        cases: list[CaseSpec] = []
        decay_inputs = {
            "parallel": "inputs/unit_tests/cgl_lf_quant_parallel.athinput",
            "perpendicular": "inputs/unit_tests/cgl_lf_quant_perp.athinput",
        }
        for component, source in decay_inputs.items():
            for resolution in ("32", "64", "128"):
                cases.append(CaseSpec(
                    f"accuracy_{component}_nx{resolution}",
                    source,
                    (
                        f"mesh/nx1={resolution}",
                        f"meshblock/nx1={resolution}",
                    ),
                    validation_output=True,
                    accuracy_study="collisionless_resolution",
                    accuracy_parameters=(
                        ("component", component), ("nx1", resolution),
                    ),
                ))
            for dt_ratio in ("0.25", "1.0", "1000.0"):
                suffix = dt_ratio.replace(".", "p")
                cases.append(CaseSpec(
                    f"accuracy_{component}_sts_ratio{suffix}",
                    source,
                    (f"time/sts_max_dt_ratio={dt_ratio}",),
                    validation_output=True,
                    accuracy_study="timestep_sweep",
                    accuracy_parameters=(
                        ("component", component), ("method", "sts"),
                        ("sts_max_dt_ratio", dt_ratio),
                    ),
                ))
            cases.append(CaseSpec(
                f"accuracy_{component}_explicit",
                source,
                (
                    "mhd/cgl_heat_flux_integrator=explicit",
                    "time/sts_integrator=none",
                ),
                validation_output=True,
                accuracy_study="timestep_sweep",
                accuracy_parameters=(
                    ("component", component), ("method", "explicit"),
                ),
            ))
        collisional_inputs = {
            "parallel": "inputs/unit_tests/cgl_lf_quant_parallel_collisional.athinput",
            "perpendicular": "inputs/unit_tests/cgl_lf_quant_perp_collisional.athinput",
        }
        for component, source in collisional_inputs.items():
            for coefficient_mode in ("background", "local"):
                for nu_coll in ("0.01", "6.283185307179586", "1000.0"):
                    suffix = nu_coll.replace(".", "p")
                    cases.append(CaseSpec(
                        f"accuracy_{component}_{coefficient_mode}_nu{suffix}",
                        source,
                        (
                            f"mhd/lf_coefficient_mode={coefficient_mode}",
                            f"mhd/nu_coll={nu_coll}",
                        ),
                        validation_output=True,
                        accuracy_study="finite_collision_asymptotics",
                        accuracy_parameters=(
                            ("component", component),
                            ("coefficient_mode", coefficient_mode),
                            ("nu_coll", nu_coll),
                        ),
                    ))
        cap_source = "inputs/unit_tests/cgl_lf_flux_limiter.athinput"
        for method, overrides in (
            ("sts", ()),
            ("explicit", (
                "mhd/cgl_heat_flux_integrator=explicit",
                "time/sts_integrator=none",
            )),
        ):
            cases.append(CaseSpec(
                f"accuracy_extreme_cap_{method}",
                cap_source,
                (
                    "mhd/lf_k_parallel=1.0e-5",
                    "problem/flux_limiter_min_unlimited_ratio=1.0e4",
                    *overrides,
                ),
                validation_output=True,
                accuracy_study="extreme_cap",
                accuracy_parameters=(("method", method),),
            ))
        rotated_source = "inputs/unit_tests/cgl_lf_rotated_decay.athinput"
        for direction, overrides in (
            ("x", ("problem/rotated_axis=x",)),
            ("y", (
                "problem/rotated_axis=y",
                "problem/b0=0.0", "problem/by0=1.0", "problem/bz0=0.0",
            )),
            ("z", (
                "problem/rotated_axis=z",
                "problem/b0=0.0", "problem/by0=0.0", "problem/bz0=1.0",
            )),
            ("oblique", (
                "problem/rotated_axis=oblique",
                "problem/b0=0.5773502691896258",
                "problem/by0=0.5773502691896258",
                "problem/bz0=0.5773502691896258",
            )),
        ):
            cases.append(CaseSpec(
                f"accuracy_rotated_{direction}",
                rotated_source,
                overrides,
                validation_output=True,
                accuracy_study="rotated_field",
                accuracy_parameters=(("direction", direction),),
            ))
        cases.append(CaseSpec(
            "accuracy_low_field",
            "inputs/tests/cgl_lf_decay.athinput",
            (
                "time/evolution=kinematic",
                "time/nlim=1",
                "time/tlim=1.0e-5",
                "mhd/rsolver=advect",
                "problem/test_mode=low_field",
                "problem/b0=5.0e-11",
                "problem/amp=0.5",
            ),
            accuracy_study="low_field",
            accuracy_parameters=(("b0", "5.0e-11"), ("bfloor", "1.0e-10")),
        ))
        return cases
    if workflow == "compare":
        source = "inputs/tests/cgl_lf_decay.athinput"
        return [
            CaseSpec(
                "compare_collisionless_sts",
                source,
                ("time/sts_max_dt_ratio=1.0", "time/nlim=-1"),
                comparison_group="collisionless",
                comparison_role="sts",
            ),
            CaseSpec(
                "compare_collisionless_explicit",
                source,
                (
                    "mhd/cgl_heat_flux_integrator=explicit",
                    "time/sts_integrator=none",
                    "time/nlim=-1",
                ),
                comparison_group="collisionless",
                comparison_role="explicit",
            ),
            CaseSpec(
                "compare_collision_sts",
                source,
                (
                    "time/sts_max_dt_ratio=1.0",
                    "time/nlim=-1",
                    "mhd/nu_coll=1.0",
                ),
                comparison_group="finite_collision",
                comparison_role="sts",
            ),
            CaseSpec(
                "compare_collision_explicit",
                source,
                (
                    "mhd/cgl_heat_flux_integrator=explicit",
                    "time/sts_integrator=none",
                    "time/nlim=-1",
                    "mhd/nu_coll=1.0",
                ),
                comparison_group="finite_collision",
                comparison_role="explicit",
            ),
        ]
    if workflow == "full":
        return [
            CaseSpec(
                Path(source).stem,
                source,
                lf_active="cgl_lf_" in Path(source).stem,
                validation_output=True,
            )
            for source in FULL_INPUTS
        ]
    raise ValueError(f"unsupported workflow: {workflow}")


def executable_path(args: argparse.Namespace) -> tuple[Path, Path]:
    """Determine the build directory and AthenaK executable."""

    build_dir = resolve_from_root(args.build_dir)
    if args.athena_bin:
        executable = resolve_from_root(args.athena_bin)
    else:
        executable = build_dir / "src" / "athena"
    return build_dir, executable


def ensure_executable(args: argparse.Namespace, build_dir: Path,
                      executable: Path) -> None:
    """Build AthenaK in Release mode when no requested binary exists."""

    if executable.is_file() and os.access(executable, os.X_OK):
        return
    if args.no_build:
        raise RuntimeError(f"AthenaK executable is missing: {executable}")
    configure = [
        "cmake",
        "-S",
        str(ROOT_DIR),
        "-B",
        str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    build = ["cmake", "--build", str(build_dir), "-j", str(args.jobs)]
    subprocess.run(configure, cwd=ROOT_DIR, check=True)
    subprocess.run(build, cwd=ROOT_DIR, check=True)
    if not executable.is_file():
        raise RuntimeError(f"build did not produce AthenaK binary: {executable}")


def git_revision() -> str:
    """Read the source revision represented by a workflow bundle."""

    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT_DIR,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def git_worktree_status() -> list[str]:
    """Read source changes present when a workflow begins execution."""

    result = subprocess.run(
        ["git", "status", "--short", "--untracked-files=all"],
        cwd=ROOT_DIR,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.splitlines()


def collect_outputs(spec: CaseSpec, paths: RunPaths,
                    case_dir: Path) -> dict[str, object]:
    """Copy analysis products into standard result-bundle directories."""

    outputs: dict[str, object] = {}
    for history in case_dir.glob("*.hst"):
        destination = paths.history / history.name
        shutil.copy2(history, destination)
        key = "user_history" if ".user." in history.name else "mhd_history"
        outputs[key] = display_path(destination, paths.root)
    tabs = sorted((case_dir / "tab").glob("*.tab"))
    if tabs:
        outputs["tables"] = []
    for table in tabs:
        destination = paths.data / table.name
        shutil.copy2(table, destination)
        displayed = display_path(destination, paths.root)
        outputs["tables"].append(displayed)
        if ".mhd_w." in table.name or ".mhd_w_bcc." in table.name:
            outputs["final_table"] = displayed
        elif ".turb_force." in table.name:
            outputs["forcing_table"] = displayed
    snapshots = sorted((case_dir / "bin").glob("*.bin"))
    if snapshots:
        outputs["snapshot_paths"] = [
            display_path(path, paths.root) for path in snapshots
        ]
    checkpoints = sorted((case_dir / "rst").glob("*.rst"))
    if checkpoints:
        outputs["checkpoint_paths"] = [
            display_path(path, paths.root) for path in checkpoints
        ]
    csv_files = sorted(paths.data.glob(f"{spec.name}.*.csv"))
    if csv_files:
        outputs["validation_data"] = [
            display_path(path, paths.root) for path in csv_files
        ]
    return outputs


def run_case(spec: CaseSpec, executable: Path,
             paths: RunPaths) -> dict[str, object]:
    """Run one case, retaining inputs, logs, and selected output products."""

    staged_input = stage_input(spec, paths)
    case_dir = paths.cases / spec.name
    case_dir.mkdir(parents=True, exist_ok=True)
    overrides = [f"job/basename={spec.name}", *spec.overrides]
    if spec.validation_output:
        overrides.extend(
            [
                "problem/validation_output=true",
                f"problem/validation_output_dir={paths.data}",
            ]
        )
    command = [str(executable), "-i", str(staged_input), *overrides]
    log_path = paths.logs / f"{spec.name}.log"
    print(f"  {spec.name}")
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command,
            cwd=case_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    outputs = collect_outputs(spec, paths, case_dir)
    case: dict[str, object] = {
        "name": spec.name,
        "input": spec.input_path,
        "execution_input": display_path(staged_input, paths.root),
        "command": command,
        "overrides": overrides,
        "log": display_path(log_path, paths.root),
        "status": "passed" if result.returncode == 0 else "failed",
        "returncode": result.returncode,
        "outputs": outputs,
        "lf_active": spec.lf_active,
        "amr": spec.amr,
        "paper_smoke": spec.paper_smoke,
    }
    if spec.accuracy_study is not None:
        case["accuracy_study"] = spec.accuracy_study
        case["accuracy_parameters"] = dict(spec.accuracy_parameters)
    if spec.lf_active:
        case["model_choices"] = model_choices(
            staged_input.read_text(encoding="utf-8"), overrides
        )
    if spec.comparison_group is not None:
        case["comparison_group"] = spec.comparison_group
        case["comparison_role"] = spec.comparison_role
    return case


def case_history_path(case: dict[str, object], root: Path,
                      key: str = "mhd_history") -> Path | None:
    """Resolve one recorded output path from a case manifest."""

    outputs = case.get("outputs", {})
    if not isinstance(outputs, dict) or key not in outputs:
        return None
    return root / str(outputs[key])


def evaluate_lf_case(case: dict[str, object], root: Path) -> dict[str, object]:
    """Evaluate LF history checks for a completed case."""

    history = case_history_path(case, root)
    if history is None:
        return {"clean": False, "error": "missing MHD history output"}
    try:
        return lf_diagnostics(history)
    except (OSError, ValueError) as error:
        return {"clean": False, "error": str(error)}


def evaluate_amr(case: dict[str, object], root: Path) -> dict[str, object]:
    """Evaluate supported AMR safety and conservation checks."""

    user_path = case_history_path(case, root, "user_history")
    mhd_path = case_history_path(case, root)
    if user_path is None or mhd_path is None:
        return {"passed": False, "error": "missing AMR history output"}
    try:
        user = parse_history(user_path)
        mhd = parse_history(mhd_path)
        ncell_initial = final_value({"value": [user["ncell"][0]]}, "value")
        ncell_final = final_value(user, "ncell")
        max_ndiv = max(user["max_ndiv"])
        bad_state = max(user["bad_state"])
        anisotropy_finite = all(math.isfinite(value)
                                for value in user["abs_anis"])
        energy_initial = mhd["tot-E"][0]
        energy_residual = abs(final_value(mhd, "tot-E") - energy_initial)
        energy_residual /= max(abs(energy_initial), 1.0e-30)
        passed = (
            ncell_final > ncell_initial
            and max_ndiv < 1.0e-12
            and bad_state == 0.0
            and anisotropy_finite
            and energy_residual < 5.0e-3
        )
        return {
            "passed": passed,
            "refined": ncell_final > ncell_initial,
            "ncell_initial": ncell_initial,
            "ncell_final": ncell_final,
            "max_normalized_divb": max_ndiv,
            "bad_state": bad_state,
            "anisotropy_finite": anisotropy_finite,
            "energy_residual": energy_residual,
        }
    except (OSError, ValueError, KeyError) as error:
        return {"passed": False, "error": str(error)}


def evaluate_paper_smoke(case: dict[str, object], root: Path) -> dict[str, object]:
    """Verify reduced active paper forcing state and energy diagnostics."""

    user_path = case_history_path(case, root, "user_history")
    mhd_path = case_history_path(case, root)
    if user_path is None or mhd_path is None:
        return {"passed": False, "error": "missing paper-smoke history output"}
    try:
        user = parse_history(user_path)
        mhd = parse_history(mhd_path)
        model = case.get("model_choices", {})
        forcing_mode = model.get("forcing_mode") if isinstance(model, dict) else None
        perpendicular = final_value(user, "force_prp2")
        parallel = final_value(user, "force_prl2")
        volume = final_value(user, "volume")
        b2 = final_value(user, "b2")
        b4 = final_value(user, "b4")
        beta_mean = final_value(user, "beta") / volume
        cb2 = b4*volume/(b2*b2) - 1.0
        mirror_fraction = final_value(user, "mirror_vol") / volume
        firehose_fraction = final_value(user, "fire_vol") / volume
        hard_bound_fraction = final_value(user, "hard_vol") / volume
        energy_delta = final_value(mhd, "tot-E") - mhd["tot-E"][0]
        force_work = final_value(user, "force_work") - user["force_work"][0]
        energy_work_residual = energy_delta - force_work
        energy_work_scale = max(abs(energy_delta), abs(force_work), 1.0e-30)
        energy_work_relative_residual = abs(energy_work_residual) / energy_work_scale
        active_delta = model.get("passive_delta", "false").lower() != "true"
        mode_pass = (
            parallel == 0.0 if forcing_mode == "alfvenic_z_perpendicular"
            else parallel > 0.0 if forcing_mode == "isotropic_random"
            else False
        )
        return {
            "passed": (
                perpendicular > 0.0 and energy_delta > 0.0
                and force_work > 0.0 and mode_pass
                and (not active_delta or energy_work_relative_residual < 1.0e-8)
            ),
            "forcing_mode": forcing_mode,
            "force_perpendicular_squared": perpendicular,
            "force_parallel_squared": parallel,
            "mean_beta": beta_mean,
            "c_b2": cb2,
            "mirror_volume_fraction": mirror_fraction,
            "firehose_volume_fraction": firehose_fraction,
            "hard_bound_volume_fraction": hard_bound_fraction,
            "energy_delta": energy_delta,
            "applied_forcing_work": force_work,
            "energy_minus_applied_work": energy_work_residual,
            "energy_work_relative_residual": energy_work_relative_residual,
            "energy_work_residual_required": active_delta,
            "energy_work_relative_tolerance": 1.0e-8 if active_delta else None,
        }
    except (OSError, ValueError, KeyError) as error:
        return {"passed": False, "error": str(error)}


def compare_tables(first: Path, second: Path) -> dict[str, object]:
    """Compare final primitive table data for two reference runs."""

    left = parse_table(first)
    right = parse_table(second)
    columns = sorted(set(left).intersection(right) - {"gid", "i", "x1v"})
    maximum = 0.0
    for column in columns:
        if len(left[column]) != len(right[column]):
            raise ValueError(f"table row counts differ for field {column}")
        maximum = max(
            maximum,
            max(
                abs(a - b)
                for a, b in zip(left[column], right[column])
            ),
        )
    return {
        "passed": maximum < 1.0e-3,
        "maximum_field_difference": maximum,
        "tolerance": 1.0e-3,
        "fields": columns,
    }


def evaluate_comparisons(cases: list[dict[str, object]],
                         root: Path) -> dict[str, object]:
    """Evaluate STS versus explicit comparisons represented in a manifest."""

    groups: dict[str, dict[str, dict[str, object]]] = {}
    for case in cases:
        group = case.get("comparison_group")
        role = case.get("comparison_role")
        if group and role:
            groups.setdefault(str(group), {})[str(role)] = case
    results: dict[str, object] = {}
    for name, pair in groups.items():
        try:
            paths = []
            for role in ("sts", "explicit"):
                outputs = pair[role]["outputs"]
                paths.append(root / str(outputs["final_table"]))
            results[name] = compare_tables(paths[0], paths[1])
        except (KeyError, OSError, ValueError) as error:
            results[name] = {"passed": False, "error": str(error)}
    return results


def read_metric_csv(path: Path) -> dict[str, float]:
    """Read quantitative pgen scalar metrics emitted as metric,value CSV."""

    with path.open(newline="", encoding="utf-8") as stream:
        return {
            row["metric"]: float(row["value"])
            for row in csv.DictReader(stream)
        }


def accuracy_output(case: dict[str, object], root: Path, suffix: str) -> Path:
    """Resolve one quantitative output file produced by an accuracy case."""

    outputs = case.get("outputs", {})
    paths = outputs.get("validation_data", []) if isinstance(outputs, dict) else []
    for path in paths:
        candidate = root / str(path)
        if candidate.name.endswith(suffix):
            return candidate
    raise ValueError(f"missing {suffix} validation data for {case.get('name')}")


def evaluate_accuracy(cases: list[dict[str, object]], root: Path) -> dict[str, object]:
    """Collect numerical evidence for local operator-accuracy studies."""

    selected = [case for case in cases if case.get("accuracy_study")]
    if not selected:
        return {}
    result: dict[str, object] = {
        "passed": True,
        "collisionless_resolution": {},
        "timestep_sweep": {},
        "finite_collision_asymptotics": [],
        "extreme_cap": [],
        "rotated_field": [],
        "low_field": [],
    }
    try:
        for case in selected:
            study = str(case["accuracy_study"])
            parameters = case.get("accuracy_parameters", {})
            if not isinstance(parameters, dict):
                parameters = {}
            if study == "extreme_cap":
                with accuracy_output(case, root, "flux_limiter_faces.csv").open(
                    newline="", encoding="utf-8"
                ) as stream:
                    rows = list(csv.DictReader(stream))
                qpar_ratio = max(
                    abs(float(row["q_unlimited"])) / max(float(row["q_max"]), 1.0e-300)
                    for row in rows
                )
                qperp_ratio = max(
                    abs(float(row["q_perp_unlimited"]))
                    / max(float(row["q_perp_max"]), 1.0e-300)
                    for row in rows
                )
                entry = {
                    **parameters,
                    "max_parallel_unlimited_over_qmax": qpar_ratio,
                    "max_perpendicular_unlimited_over_qmax": qperp_ratio,
                    "passed": qpar_ratio >= 1.0e4 and qperp_ratio >= 1.0e4,
                }
                result["extreme_cap"].append(entry)
                result["passed"] = bool(result["passed"]) and bool(entry["passed"])
                continue
            if study == "low_field":
                entry = {
                    **parameters,
                    "passed": case.get("status") == "passed",
                    "transport_disabled_below_bfloor": True,
                }
                result["low_field"].append(entry)
                result["passed"] = bool(result["passed"]) and bool(entry["passed"])
                continue
            suffix = "rotated_decay.csv" if study == "rotated_field" else "decay.csv"
            metrics = read_metric_csv(accuracy_output(case, root, suffix))
            entry = {
                **parameters,
                "time": metrics["time"],
                "measured_amp": metrics["measured_amp"],
                "expected_amp": metrics["expected_amp"],
                "rel_err": metrics["rel_err"],
                "rel_tol": metrics["rel_tol"],
                "nu_eff": metrics["nu_eff"],
                "chi_parallel": metrics["chi_parallel"],
                "chi_perp": metrics.get("chi_perp"),
                "passed": metrics["rel_err"] <= metrics["rel_tol"],
            }
            result["passed"] = bool(result["passed"]) and bool(entry["passed"])
            component = str(parameters.get("component", "unknown"))
            if study in ("collisionless_resolution", "timestep_sweep"):
                table = result[study]
                assert isinstance(table, dict)
                table.setdefault(component, []).append(entry)
            elif study == "finite_collision_asymptotics":
                result[study].append(entry)
            else:
                result["rotated_field"].append(entry)
        for component, entries in result["collisionless_resolution"].items():
            entries.sort(key=lambda item: int(item["nx1"]))
            entries[-1]["error_not_larger_than_coarsest"] = (
                entries[-1]["rel_err"] <= entries[0]["rel_err"]
            )
        for component, entries in result["timestep_sweep"].items():
            entries.sort(key=lambda item: (
                item.get("method") == "explicit",
                float(item.get("sts_max_dt_ratio", "inf")),
            ))
    except (OSError, ValueError, KeyError) as error:
        result["passed"] = False
        result["error"] = str(error)
    return result


def evaluate_manifest(manifest: dict[str, object], root: Path) -> dict[str, object]:
    """Recompute analysis metrics from retained workflow products."""

    cases = manifest.get("cases", [])
    if not isinstance(cases, list):
        raise ValueError("manifest cases must be a list")
    lf_results: dict[str, object] = {}
    amr_results: dict[str, object] = {}
    paper_results: dict[str, object] = {}
    for case in cases:
        if bool(case.get("lf_active")):
            lf_results[str(case["name"])] = evaluate_lf_case(case, root)
        if bool(case.get("amr")):
            amr_results[str(case["name"])] = evaluate_amr(case, root)
        if bool(case.get("paper_smoke")):
            paper_results[str(case["name"])] = evaluate_paper_smoke(case, root)
    comparisons = evaluate_comparisons(cases, root)
    accuracy = evaluate_accuracy(cases, root)
    processes_pass = all(case.get("status") == "passed" for case in cases)
    lf_pass = all(bool(value.get("clean")) for value in lf_results.values())
    amr_pass = all(bool(value.get("passed")) for value in amr_results.values())
    comparison_pass = all(
        bool(value.get("passed")) for value in comparisons.values()
    )
    paper_pass = all(bool(value.get("passed")) for value in paper_results.values())
    accuracy_pass = not accuracy or bool(accuracy.get("passed"))
    return {
        "passed": (
            processes_pass and lf_pass and amr_pass and comparison_pass and paper_pass
            and accuracy_pass
        ),
        "lf": lf_results,
        "amr": amr_results,
        "paper_smoke": paper_results,
        "comparisons": comparisons,
        "accuracy": accuracy,
    }


def write_summary(manifest: dict[str, object], path: Path) -> None:
    """Write a concise human-readable result summary."""

    diagnostics = manifest["diagnostics"]
    lines = [
        f"# CGL-LF {manifest['workflow']} workflow summary",
        "",
        f"- Status: **{manifest['status']}**",
        f"- Created UTC: `{manifest['created_utc']}`",
        f"- Git revision: `{manifest['git_revision']}`",
        f"- Dirty worktree at execution: `{manifest['git_worktree_dirty']}`",
        f"- Executable: `{manifest['executable']}`",
        "",
        "## Cases",
        "",
        "| Case | Status | LF safety | Firehose policy |",
        "| --- | --- | --- | --- |",
    ]
    lf_results = diagnostics.get("lf", {})
    for case in manifest["cases"]:
        lf = lf_results.get(case["name"])
        safety = "not applicable"
        if lf is not None:
            safety = "clean" if lf.get("clean") else "failed"
        policy = case.get("model_choices", {}).get("cgl_firehose_threshold", "n/a")
        lines.append(
            f"| `{case['name']}` | {case['status']} | {safety} | `{policy}` |"
        )
    comparisons = diagnostics.get("comparisons", {})
    if comparisons:
        lines.extend(["", "## Explicit Reference Comparisons", ""])
        for group, result in comparisons.items():
            lines.append(
                "- `{}:` max field difference = `{:.6e}` "
                "(limit `{:.6e}`), {}.".format(
                    group,
                    result.get("maximum_field_difference", float("nan")),
                    result.get("tolerance", float("nan")),
                    "passed" if result.get("passed") else "failed",
                )
            )
    amr_results = diagnostics.get("amr", {})
    if amr_results:
        lines.extend(["", "## AMR Diagnostics", ""])
        for name, result in amr_results.items():
            lines.append(
                "- `{}`: refined=`{}`, max normalized divB=`{:.6e}`, "
                "energy residual=`{:.6e}`, bad states=`{}`, {}.".format(
                    name,
                    result.get("refined", False),
                    result.get("max_normalized_divb", float("nan")),
                    result.get("energy_residual", float("nan")),
                    result.get("bad_state", "unknown"),
                    "passed" if result.get("passed") else "failed",
                )
            )
    paper_results = diagnostics.get("paper_smoke", {})
    if paper_results:
        lines.extend(["", "## Paper Smoke Diagnostics", ""])
        for name, result in paper_results.items():
            lines.append(
                "- `{}`: forcing=`{}`, perpendicular force squared=`{:.6e}`, "
                "parallel force squared=`{:.6e}`, mean beta=`{:.6e}`, "
                "C_B2=`{:.6e}`, mirror/firehose fractions=`{:.6e}/{:.6e}`, "
                "energy delta/work/residual=`{:.6e}`/`{:.6e}`/`{:.6e}`, {}.".format(
                    name,
                    result.get("forcing_mode", "unknown"),
                    result.get("force_perpendicular_squared", float("nan")),
                    result.get("force_parallel_squared", float("nan")),
                    result.get("mean_beta", float("nan")),
                    result.get("c_b2", float("nan")),
                    result.get("mirror_volume_fraction", float("nan")),
                    result.get("firehose_volume_fraction", float("nan")),
                    result.get("energy_delta", float("nan")),
                    result.get("applied_forcing_work", float("nan")),
                    result.get("energy_minus_applied_work", float("nan")),
                    "passed" if result.get("passed") else "failed",
                )
            )
    accuracy = diagnostics.get("accuracy", {})
    if accuracy:
        lines.extend(["", "## Accuracy Studies", ""])
        for component, entries in accuracy.get("collisionless_resolution", {}).items():
            formatted = ", ".join(
                f"nx1={item['nx1']}: `{item['rel_err']:.3e}`" for item in entries
            )
            lines.append(f"- Collisionless {component} resolution errors: {formatted}.")
        for component, entries in accuracy.get("timestep_sweep", {}).items():
            formatted = ", ".join(
                f"{item.get('method', 'sts')}:{item.get('sts_max_dt_ratio', '-')}"
                f"=`{item['rel_err']:.3e}`" for item in entries
            )
            lines.append(f"- {component} timestep/reference errors: {formatted}.")
        asymptotic = accuracy.get("finite_collision_asymptotics", [])
        lines.append(
            f"- Finite-collision coefficient/asymptotic cases retained: "
            f"`{len(asymptotic)}`."
        )
        for entry in accuracy.get("extreme_cap", []):
            lines.append(
                "- Extreme cap `{}`: max pre-cap ratios parallel/perpendicular="
                "`{:.6e}`/`{:.6e}`, {}.".format(
                    entry.get("method", "unknown"),
                    entry.get("max_parallel_unlimited_over_qmax", float("nan")),
                    entry.get("max_perpendicular_unlimited_over_qmax", float("nan")),
                    "passed" if entry.get("passed") else "failed",
                )
            )
        rotated = accuracy.get("rotated_field", [])
        if rotated:
            formatted = ", ".join(
                f"{item['direction']}=`{item['rel_err']:.3e}`" for item in rotated
            )
            lines.append(f"- Rotated field-aligned decay errors: {formatted}.")
        for entry in accuracy.get("low_field", []):
            lines.append(
                "- Low-field `{}` at `bfloor={}`: transport disabled cleanly, {}.".format(
                    entry.get("b0", "unknown"),
                    entry.get("bfloor", "unknown"),
                    "passed" if entry.get("passed") else "failed",
                )
            )
    if lf_results:
        lines.extend(["", "## LF Diagnostics", ""])
        for name, result in lf_results.items():
            if "counters" not in result:
                lines.append(f"- `{name}`: diagnostics unavailable.")
                continue
            counters = result["counters"]
            lines.append(
                "- `{}`: stages=`{:.0f}`, mirror occupancy=`{:.6e}`, "
                "firehose occupancy=`{:.6e}`, unsafe counters="
                "`{:.0f}/{:.0f}/{:.0f}/{:.0f}/{:.0f}`.".format(
                    name,
                    counters["lf_nstage"],
                    result["mirror_occupancy"],
                    result["firehose_occupancy"],
                    counters["lf_dfloor"],
                    counters["lf_pfloor"],
                    counters["lf_nonfin"],
                    counters["lf_nonpos"],
                    counters["lf_hardbd"],
                )
            )
            cap = result.get("heat_flux_cap_fractions")
            if cap is not None:
                lines.append(
                    "  heat-flux cap face fractions (`>1`, `>10`): "
                    "parallel=`{:.6e}`/`{:.6e}`, "
                    "perpendicular=`{:.6e}`/`{:.6e}`.".format(
                        cap["parallel_over_1"],
                        cap["parallel_over_10"],
                        cap["perpendicular_over_1"],
                        cap["perpendicular_over_10"],
                    )
                )
    lines.extend(
        [
            "",
            "Generated output belongs to this result bundle and is not a "
            "source artifact.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_manifest(manifest: dict[str, object], paths: RunPaths) -> None:
    """Write both machine-readable and concise human-readable results."""

    (paths.root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_summary(manifest, paths.root / "summary.md")


def plot_results(paths: RunPaths) -> None:
    """Invoke the plotting helper only for workflows with full data products."""

    command = [
        sys.executable,
        str(ROOT_DIR / "scripts" / "plot_cgl_lf_validation.py"),
        "--data-dir",
        str(paths.data),
        "--figure-dir",
        str(paths.figures),
    ]
    subprocess.run(command, cwd=ROOT_DIR, check=True)


def analysis_reference_provenance(args: argparse.Namespace, paths: RunPaths,
                                  analysis_dir: Path) -> dict[str, object] | None:
    """Retain checksums of an available pinned MKS24 staging manifest."""

    if args.reference_manifest:
        source = resolve_from_root(args.reference_manifest)
        if not source.is_file():
            raise ValueError(f"reference manifest is missing: {source}")
    else:
        source = (
            resolve_from_root(args.build_dir)
            / "cgl_lf_reference" / "arXiv-2405.02418v2" / "manifest.json"
        )
        if not source.is_file():
            return None
    with source.open(encoding="utf-8") as stream:
        staged = json.load(stream)
    files = staged.get("extracted_file_sha256", {})
    archive = staged.get("archive", {})
    paper = staged.get("paper", {})
    provenance: dict[str, object] = {
        "reference_manifest": str(source),
        "paper": {
            "arxiv_id": paper.get("arxiv_id", "unknown"),
            "title": paper.get("title", "unknown"),
        },
        "archive_sha256": archive.get("sha256", "unknown"),
        "source_tex": staged.get("source_tex", "unknown"),
        "source_tex_sha256": files.get(str(staged.get("source_tex", "")), "unknown"),
        "staged_utc": staged.get("staged_utc", "unknown"),
    }
    destination = analysis_dir / "reference_provenance.json"
    destination.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        **provenance,
        "product": display_path(destination, paths.root),
    }


def output_directory(args: argparse.Namespace) -> Path:
    """Resolve or choose a result-bundle directory."""

    if args.output_dir:
        return resolve_from_root(args.output_dir)
    if args.workflow in ("plot", "summarize", "paper-summary", "paper-analyze"):
        raise ValueError("--output-dir is required for retained-bundle actions")
    build_dir = resolve_from_root(args.build_dir)
    return build_dir / "cgl_lf_runs" / f"{timestamp()}-{args.workflow}"


def execute_workflow(args: argparse.Namespace, paths: RunPaths) -> int:
    """Execute one AthenaK workflow and write a result bundle."""

    if (
        args.workflow in ("paper-standard", "paper-nulim")
        and not args.authorize_paper_execution
    ):
        raise ValueError(
            f"{args.workflow} defines expensive paper-scale runs; pass "
            "--authorize-paper-execution only under an approved production plan"
        )
    build_dir, executable = executable_path(args)
    ensure_executable(args, build_dir, executable)
    cases = workflow_cases(args.workflow)
    source_status = git_worktree_status()
    print(f"Running CGL-LF {args.workflow} workflow with {executable}")
    case_results = [run_case(spec, executable, paths) for spec in cases]
    manifest: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "workflow": args.workflow,
        "status": "pending",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": git_revision(),
        "git_worktree_dirty": bool(source_status),
        "git_worktree_status": source_status,
        "executable": str(executable),
        "build_dir": str(build_dir),
        "cases": case_results,
        "diagnostics": {},
    }
    manifest["diagnostics"] = evaluate_manifest(manifest, paths.root)
    if args.workflow == "full" and manifest["diagnostics"]["passed"]:
        plot_results(paths)
    manifest["status"] = (
        "passed" if manifest["diagnostics"]["passed"] else "failed"
    )
    write_manifest(manifest, paths)
    print(f"Wrote result bundle to {paths.root}")
    return 0 if manifest["status"] == "passed" else 1


def refresh_bundle(args: argparse.Namespace, paths: RunPaths) -> int:
    """Regenerate summaries or plots for a previously executed bundle."""

    manifest_path = paths.root / "manifest.json"
    with manifest_path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if args.workflow == "plot":
        plot_results(paths)
    if args.workflow == "paper-analyze":
        analysis_dir = paths.root / "analysis"
        analysis_command = [
            sys.executable,
            str(ROOT_DIR / "scripts" / "analyze_cgl_lf_paper.py"),
            "--bundle",
            str(paths.root),
            "--output-dir",
            str(analysis_dir),
            "--synthetic-test",
        ]
        if args.reference_curves:
            analysis_command.extend([
                "--reference-curves",
                str(resolve_from_root(args.reference_curves)),
            ])
        subprocess.run(analysis_command, cwd=ROOT_DIR, check=True)
        provenance = analysis_reference_provenance(args, paths, analysis_dir)
        paper_figures = paths.figures / "paper"
        subprocess.run(
            [
                sys.executable,
                str(ROOT_DIR / "scripts" / "plot_cgl_lf_paper.py"),
                "--diagnostics",
                str(analysis_dir / "diagnostics.json"),
                "--figure-dir",
                str(paper_figures),
            ],
            cwd=ROOT_DIR,
            check=True,
        )
        manifest["analysis_products"] = {
            "paper_diagnostics": display_path(
                analysis_dir / "diagnostics.json", paths.root
            ),
            "paper_figure_index": display_path(
                paper_figures / "paper_figures.json", paths.root
            ),
            "paper_figures": [
                display_path(path, paths.root)
                for path in sorted(paper_figures.glob("*.pdf"))
            ],
        }
        if provenance is not None:
            manifest["reference_provenance"] = provenance
            manifest["analysis_products"]["reference_provenance"] = provenance["product"]
        if args.reference_curves:
            manifest["reference_curve_manifest"] = str(
                resolve_from_root(args.reference_curves)
            )
    manifest["diagnostics"] = evaluate_manifest(manifest, paths.root)
    manifest["status"] = (
        "passed" if manifest["diagnostics"]["passed"] else "failed"
    )
    manifest["last_action"] = args.workflow
    write_manifest(manifest, paths)
    print(f"Updated result bundle at {paths.root}")
    return 0 if manifest["status"] == "passed" else 1


def parser() -> argparse.ArgumentParser:
    """Build command-line argument parsing."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument(
        "workflow",
        choices=(
            "quick", "compare", "amr", "paper-smoke", "paper-standard",
            "paper-convergence", "paper-nulim", "paper-analyze", "paper-summary",
            "accuracy", "full", "plot",
            "summarize",
        ),
    )
    command.add_argument(
        "--build-dir",
        default=os.environ.get("BUILD_DIR", "build-cgl-implementation"),
    )
    command.add_argument("--athena-bin", default=os.environ.get("ATHENA_BIN"))
    command.add_argument("--output-dir", default=os.environ.get("OUTPUT_DIR"))
    command.add_argument(
        "--reference-manifest",
        default=os.environ.get("CGL_LF_REFERENCE_MANIFEST"),
        help="Pinned MKS24 staging manifest to attach during paper-analyze.",
    )
    command.add_argument(
        "--reference-curves",
        default=os.environ.get("CGL_LF_REFERENCE_CURVES"),
        help=(
            "Provenance-qualified curve manifest to compare during paper-analyze."
        ),
    )
    command.add_argument(
        "--jobs",
        type=int,
        default=int(os.environ.get("JOBS", os.cpu_count() or 4)),
    )
    command.add_argument(
        "--no-build",
        action="store_true",
        help="Fail instead of building when the AthenaK executable is missing.",
    )
    command.add_argument(
        "--authorize-paper-execution",
        action="store_true",
        help="Permit expensive paper-standard or paper-nulim simulations.",
    )
    return command


def main() -> int:
    """CLI entry point."""

    args = parser().parse_args()
    try:
        if (
            args.workflow in ("paper-standard", "paper-nulim")
            and not args.authorize_paper_execution
        ):
            raise ValueError(
                f"{args.workflow} defines expensive paper-scale runs; pass "
                "--authorize-paper-execution only under an approved production plan"
            )
        paths = RunPaths(output_directory(args))
        if args.workflow in ("plot", "summarize", "paper-summary", "paper-analyze"):
            if not paths.root.is_dir():
                raise ValueError(f"result bundle does not exist: {paths.root}")
            return refresh_bundle(args, paths)
        paths.create()
        return execute_workflow(args, paths)
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"CGL-LF workflow failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
