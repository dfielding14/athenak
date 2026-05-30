#!/usr/bin/env python3
"""Exercise retained relativistic diagnostic writers at numerical boundaries."""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Optional, Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_INPUT = REPO_ROOT / "inputs" / "particles" / (
    "cr_tracer_relativistic_contract.athinput"
)
sys.path.insert(0, str(SCRIPT_DIR))

import cr_tracer_inspect as inspector  # noqa: E402


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _pass(label: str, detail: str = "") -> None:
    suffix = f": {detail}" if detail else ""
    print(f"PASS: {label}{suffix}")


def _flags(**values: object) -> list[str]:
    return [f"{key.replace('__', '/')}={value}" for key, value in values.items()]


def _output_block(file_type: str, extra: str = "") -> str:
    return (
        "\n<output1>\n"
        f"file_type = {file_type}\n"
        "dcycle    = 1\n"
        f"{extra}"
    )


def _write_input(case_dir: Path, base_input: Path, file_type: str,
                 extra: str = "") -> Path:
    input_path = case_dir / "case.athinput"
    input_path.write_text(base_input.read_text() + _output_block(file_type, extra))
    return input_path


def _run(binary: Path, case_dir: Path, input_path: Path,
         flags: Sequence[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [str(binary), "-i", str(input_path), *flags],
        cwd=case_dir,
        capture_output=True,
        text=True,
        timeout=60.0,
    )
    (case_dir / "athena.log").write_text(result.stdout + result.stderr)
    return result


def _expect_failure(result: subprocess.CompletedProcess[str],
                    substring: str, label: str) -> None:
    output = result.stdout + result.stderr
    _require(result.returncode != 0, f"{label}: Athena unexpectedly succeeded")
    _require(substring in output, f"{label}: missing failure text {substring!r}")


def _common_flags(basename: str) -> list[str]:
    return _flags(
        job__basename=basename,
        problem__particle_position="fixed",
        problem__particle_x1=0.0,
        problem__particle_x2=0.0,
        problem__particle_x3=0.0,
        problem__particle_velocity="uniform",
        problem__v0x=0.2,
        problem__v0y=0.0,
        problem__v0z=0.0,
    )


def _pmom_scaled_norm_success(binary: Path, root: Path, base_input: Path) -> None:
    case_dir = root / "pmom_scaled_norm_success"
    case_dir.mkdir()
    input_path = _write_input(case_dir, base_input, "pmom")
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("pmom_scaled_norm_success") + _flags(
            problem__B0x=1.0e308,
            problem__B0y=0.0,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "pmom scaled norm: Athena failed")
    record = inspector.read_pmom_record(
        case_dir / "pmom" / "pmom_scaled_norm_success.pmom",
        nspecies=1,
    )
    mean_mu = record["moments"][0]["mean_mu"]
    _require(
        math.isclose(mean_mu, 1.0, rel_tol=0.0, abs_tol=1.0e-12),
        f"pmom scaled norm: expected mean_mu=1, got {mean_mu}",
    )


def _pmom_preaccumulation_rejection(binary: Path, root: Path,
                                    base_input: Path) -> None:
    case_dir = root / "pmom_preaccumulation_rejection"
    case_dir.mkdir()
    input_path = _write_input(case_dir, base_input, "pmom")
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("pmom_preaccumulation_rejection") + _flags(
            particles__c_model=1.0e155,
            problem__v0x=1.0e154,
            problem__v0y=1.0e154,
            problem__v0z=1.0e154,
        ),
    )
    _expect_failure(
        result,
        "relativistic_hc pmom derived quantity evaluation failed",
        "pmom pre-accumulation",
    )
    _require(
        not list((case_dir / "pmom").glob("*.pmom")),
        "pmom pre-accumulation: writer left a partial artifact",
    )


def _df_scaled_norm_success(binary: Path, root: Path, base_input: Path) -> None:
    case_dir = root / "df_scaled_norm_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "df",
        "nbin = 10\n"
        "vmin = 0.05\n"
        "vmax = 0.15\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("df_scaled_norm_success") + _flags(
            particles__c_model=10.0,
            problem__v0x=1.0,
            problem__v0y=-0.8,
            problem__v0z=0.0,
            problem__B0x=1.0e308,
            problem__B0y=1.0e308,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "df scaled norm: Athena failed")
    histogram = inspector.read_df_file(
        case_dir / "df" / "df_scaled_norm_success.df",
        nspecies=1,
        nbin=10,
    )
    _require(
        histogram == [[0, 0, 0, 0, 0, 0, 1, 0, 0, 0]],
        f"df scaled norm: unexpected histogram {histogram}",
    )


def _pspec_mu_scaled_speed_success(binary: Path, root: Path,
                                   base_input: Path) -> None:
    case_dir = root / "pspec_mu_scaled_speed_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "pspec",
        "quantity = mu\n"
        "nbin = 10\n"
        "vmin = 0.65\n"
        "vmax = 0.75\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("pspec_mu_scaled_speed_success") + _flags(
            particles__c_model=1.0e155,
            problem__v0x=1.0e154,
            problem__v0y=1.0e154,
            problem__v0z=0.0,
            problem__B0x=1.0,
            problem__B0y=0.0,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "pspec scaled speed: Athena failed")
    histogram = inspector.read_pspec_file(
        case_dir / "pspec" / "pspec_mu_scaled_speed_success.pspec",
        nspecies=1,
        nbin=10,
    )
    _require(
        histogram == [[0, 0, 0, 0, 0, 1, 0, 0, 0, 0]],
        f"pspec scaled speed: unexpected histogram {histogram}",
    )


def _pspec_proxy_stability_success(binary: Path, root: Path,
                                   base_input: Path) -> None:
    case_dir = root / "pspec_proxy_stability_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "pspec",
        "quantity = velocity_magnetic_moment_proxy\n"
        "nbin = 10\n"
        "vmin = 0.0\n"
        "vmax = 2e-24\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("pspec_proxy_stability_success") + _flags(
            problem__v0x=0.2,
            problem__v0y=1.0e-12,
            problem__v0z=0.0,
            problem__B0x=1.0,
            problem__B0y=0.0,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "pspec stable proxy: Athena failed")
    histogram = inspector.read_pspec_file(
        case_dir / "pspec" / "pspec_proxy_stability_success.pspec",
        nspecies=1,
        nbin=10,
    )
    _require(
        histogram == [[0, 0, 0, 0, 0, 1, 0, 0, 0, 0]],
        f"pspec stable proxy: unexpected histogram {histogram}",
    )


def _pspec_proxy_tiny_b_success(binary: Path, root: Path,
                                base_input: Path) -> None:
    case_dir = root / "pspec_proxy_tiny_b_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "pspec",
        "quantity = velocity_magnetic_moment_proxy\n"
        "nbin = 10\n"
        "vmin = 1e303\n"
        "vmax = 1.02e303\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("pspec_proxy_tiny_b_success") + _flags(
            problem__v0x=0.0,
            problem__v0y=1.0e-10,
            problem__v0z=0.0,
            problem__B0x=1.0e-323,
            problem__B0y=0.0,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "pspec tiny-B proxy: Athena failed")
    histogram = inspector.read_pspec_file(
        case_dir / "pspec" / "pspec_proxy_tiny_b_success.pspec",
        nspecies=1,
        nbin=10,
    )
    _require(
        histogram == [[0, 0, 0, 0, 0, 0, 1, 0, 0, 0]],
        f"pspec tiny-B proxy: unexpected histogram {histogram}",
    )


def _pspec2_scaled_speed_success(binary: Path, root: Path,
                                 base_input: Path) -> None:
    case_dir = root / "pspec2_scaled_speed_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "pspec2",
        "quantity = mu_wmag\n"
        "nbin1 = 10\n"
        "vmin1 = 0.65\n"
        "vmax1 = 0.75\n"
        "nbin2 = 10\n"
        "vmin2 = 1e154\n"
        "vmax2 = 2e154\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("pspec2_scaled_speed_success") + _flags(
            particles__c_model=1.0e155,
            problem__v0x=1.0e154,
            problem__v0y=1.0e154,
            problem__v0z=0.0,
            problem__B0x=1.0,
            problem__B0y=0.0,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "pspec2 scaled speed: Athena failed")
    histogram = inspector.read_pspec2_file(
        case_dir / "pspec2" / "pspec2_scaled_speed_success.pspec2",
        nspecies=1,
        nbin1=10,
        nbin2=10,
    )
    _require(
        histogram[0][5][4] == 1,
        f"pspec2 scaled speed: unexpected histogram {histogram}",
    )


def _pspec2_vperp_stability_success(binary: Path, root: Path,
                                    base_input: Path) -> None:
    case_dir = root / "pspec2_vperp_stability_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "pspec2",
        "quantity = vpar_vperp\n"
        "nbin1 = 10\n"
        "vmin1 = 0.15\n"
        "vmax1 = 0.25\n"
        "nbin2 = 10\n"
        "vmin2 = 0.0\n"
        "vmax2 = 2e-12\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("pspec2_vperp_stability_success") + _flags(
            problem__v0x=0.2,
            problem__v0y=1.0e-12,
            problem__v0z=0.0,
            problem__B0x=1.0,
            problem__B0y=0.0,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "pspec2 stable vperp: Athena failed")
    histogram = inspector.read_pspec2_file(
        case_dir / "pspec2" / "pspec2_vperp_stability_success.pspec2",
        nspecies=1,
        nbin1=10,
        nbin2=10,
    )
    _require(
        histogram[0][5][5] == 1,
        f"pspec2 stable vperp: unexpected histogram {histogram}",
    )


def _psamp_vperp_stability_success(binary: Path, root: Path,
                                   base_input: Path) -> None:
    case_dir = root / "psamp_vperp_stability_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "psamp",
        "species = -1\n"
        "sample_stride = 1\n"
        "sample_offset = 0\n"
        "fields = vperp,velocity_magnetic_moment_proxy\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("psamp_vperp_stability_success") + _flags(
            problem__v0x=0.2,
            problem__v0y=1.0e-12,
            problem__v0z=0.0,
            problem__B0x=1.0,
            problem__B0y=0.0,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "psamp stable vperp: Athena failed")
    sample = inspector.read_psamp_file(
        case_dir / "psamp" / "rank_00000000" /
        "psamp_vperp_stability_success.psamp"
    )
    _require(len(sample["rows"]) == 1, "psamp stable vperp: expected one row")
    fields = sample["rows"][0]["fields"]
    _require(
        math.isclose(fields["vperp"], 1.0e-12, rel_tol=1.0e-12, abs_tol=0.0),
        f"psamp stable vperp: unexpected vperp {fields['vperp']}",
    )
    _require(
        math.isclose(
            fields["velocity_magnetic_moment_proxy"],
            1.0e-24,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ),
        "psamp stable vperp: unexpected magnetic-moment proxy "
        f"{fields['velocity_magnetic_moment_proxy']}",
    )


def _psamp_proxy_tiny_b_success(binary: Path, root: Path,
                                base_input: Path) -> None:
    case_dir = root / "psamp_proxy_tiny_b_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "psamp",
        "species = -1\n"
        "sample_stride = 1\n"
        "sample_offset = 0\n"
        "fields = velocity_magnetic_moment_proxy\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("psamp_proxy_tiny_b_success") + _flags(
            problem__v0x=0.0,
            problem__v0y=1.0e-10,
            problem__v0z=0.0,
            problem__B0x=1.0e-323,
            problem__B0y=0.0,
            problem__B0z=0.0,
        ),
    )
    _require(result.returncode == 0, "psamp tiny-B proxy: Athena failed")
    sample = inspector.read_psamp_file(
        case_dir / "psamp" / "rank_00000000" /
        "psamp_proxy_tiny_b_success.psamp"
    )
    _require(len(sample["rows"]) == 1, "psamp tiny-B proxy: expected one row")
    proxy = sample["rows"][0]["fields"]["velocity_magnetic_moment_proxy"]
    expected = (1.0e-10 * 1.0e-10) / float("1e-323")
    _require(
        math.isclose(proxy, expected, rel_tol=1.0e-12, abs_tol=0.0),
        f"psamp tiny-B proxy: expected {expected}, got {proxy}",
    )


def _psamp_row_atomic_rejection(binary: Path, root: Path,
                                base_input: Path) -> None:
    case_dir = root / "psamp_row_atomic_rejection"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "psamp",
        "species = -1\n"
        "sample_stride = 1\n"
        "sample_offset = 0\n"
        "fields = x,r_larmor_over_dx_min\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("psamp_row_atomic_rejection") + _flags(
            problem__B0x=0.0,
            problem__B0y=0.0,
            problem__B0z=5.0e-310,
        ),
    )
    _expect_failure(
        result,
        "relativistic_hc psamp derived field evaluation failed",
        "psamp row atomicity",
    )
    psamp_path = (
        case_dir / "psamp" / "rank_00000000" /
        "psamp_row_atomic_rejection.psamp"
    )
    _require(psamp_path.is_file(), "psamp row atomicity: missing header artifact")
    data_rows = [
        line for line in psamp_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    _require(not data_rows, "psamp row atomicity: writer left a partial data row")


def _psamp_ratio_scaled_success(binary: Path, root: Path,
                                base_input: Path) -> None:
    case_dir = root / "psamp_ratio_scaled_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "psamp",
        "species = -1\n"
        "sample_stride = 1\n"
        "sample_offset = 0\n"
        "fields = r_larmor_over_dx_min\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("psamp_ratio_scaled_success") + _flags(
            mesh__x1min=-5.0e307,
            mesh__x1max=5.0e307,
            mesh__x2min=-5.0e307,
            mesh__x2max=5.0e307,
            mesh__x3min=-5.0e307,
            mesh__x3max=5.0e307,
            particles__alpha_s=20.0,
        ),
    )
    _require(result.returncode == 0, "psamp scaled ratio: Athena failed")
    psamp_path = (
        case_dir / "psamp" / "rank_00000000" /
        "psamp_ratio_scaled_success.psamp"
    )
    data_rows = [
        line for line in psamp_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    _require(data_rows, "psamp scaled ratio: expected at least one data row")
    gamma = 1.0 / math.sqrt(1.0 - 0.2 * 0.2)
    expected = gamma * 0.2 / 20.0 / 1.25e307
    for row in data_rows:
        ratio = float(row.split()[-1])
        _require(
            math.isclose(ratio, expected, rel_tol=1.0e-12, abs_tol=0.0),
            f"psamp scaled ratio: expected {expected}, got {ratio}",
        )


def _zero_field_rejections(binary: Path, root: Path, base_input: Path) -> None:
    cases = (
        ("df", "", "particle df derived quantity evaluation failed"),
        ("pmom", "", "relativistic_hc pmom derived quantity evaluation failed"),
        (
            "pspec",
            "quantity = mu\nnbin = 8\nvmin = -1.0\nvmax = 1.0\n",
            "relativistic_hc pspec derived quantity evaluation failed",
        ),
        (
            "pspec2",
            "quantity = mu_wmag\n"
            "nbin1 = 8\nvmin1 = -1.0\nvmax1 = 1.0\n"
            "nbin2 = 8\nvmin2 = 0.0\nvmax2 = 1.0\n",
            "relativistic_hc pspec2 derived quantity evaluation failed",
        ),
        (
            "psamp",
            "species = -1\nsample_stride = 1\nsample_offset = 0\n"
            "fields = r_larmor_over_dx_min\n",
            "relativistic_hc psamp derived field evaluation failed",
        ),
    )
    for file_type, extra, failure in cases:
        case_dir = root / f"{file_type}_zero_field_rejection"
        case_dir.mkdir()
        input_path = _write_input(case_dir, base_input, file_type, extra)
        result = _run(
            binary,
            case_dir,
            input_path,
            _common_flags(f"{file_type}_zero_field_rejection") + _flags(
                problem__B0x=0.0,
                problem__B0y=0.0,
                problem__B0z=0.0,
            ),
        )
        _expect_failure(result, failure, f"{file_type} zero field")
        _pass(f"{file_type} zero-field pre-artifact rejection")


def _ppd_preopen_rejection(binary: Path, root: Path, base_input: Path) -> None:
    case_dir = root / "ppd_preopen_rejection"
    case_dir.mkdir()
    input_path = _write_input(case_dir, base_input, "ppd")
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("ppd_preopen_rejection") + _flags(
            mesh__x1min=-2.0e39,
            mesh__x1max=2.0e39,
            problem__particle_x1=1.5e39,
        ),
    )
    _expect_failure(
        result,
        "relativistic_hc ppd values must be finite and identities exactly "
        "representable in float32",
        "ppd pre-open",
    )
    _require(
        not list((case_dir / "ppd").glob("*.ppd")),
        "ppd pre-open: writer left a partial artifact",
    )
    _pass("ppd non-finite float32 position pre-artifact rejection")


def _pspec_log_floor_metadata_success(binary: Path, root: Path,
                                      base_input: Path) -> None:
    case_dir = root / "pspec_log_floor_metadata_success"
    case_dir.mkdir()
    input_path = _write_input(
        case_dir,
        base_input,
        "pspec",
        "quantity = log10_kinetic_energy_model\n"
        "nbin = 10\n"
        "vmin = -400\n"
        "vmax = 1\n",
    )
    result = _run(
        binary,
        case_dir,
        input_path,
        _common_flags("pspec_log_floor_metadata_success") + _flags(
            problem__v0x=0.0,
            problem__v0y=0.0,
            problem__v0z=0.0,
        ),
    )
    _require(result.returncode == 0, "pspec logarithmic floor metadata: Athena failed")
    record = inspector.read_pspec_record(
        case_dir / "pspec" / "pspec_log_floor_metadata_success.pspec",
        nspecies=1,
        nbin=10,
    )
    metadata = record["metadata"]
    _require(
        metadata.get("log10_floor_units") == "code_velocity_squared",
        "pspec logarithmic floor metadata: missing declared floor units",
    )
    _require(
        math.isclose(float(metadata.get("log10_floor", "nan")), 1.0e-30,
                     rel_tol=0.0, abs_tol=0.0),
        f"pspec logarithmic floor metadata: unexpected floor "
        f"{metadata.get('log10_floor')!r}",
    )
    _require(
        record["histogram"][0] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        f"pspec logarithmic floor metadata: unexpected histogram "
        f"{record['histogram']}",
    )
    _pass(
        "pspec zero-state logarithmic-floor writer artifact",
        "log10_floor=1e-30 log10_floor_units=code_velocity_squared occupied_bin=9",
    )


@contextmanager
def _runtime_root(work_dir: Optional[Path]) -> Iterator[Path]:
    if work_dir is not None:
        work_dir.mkdir(parents=True, exist_ok=True)
        yield work_dir
        return
    with tempfile.TemporaryDirectory(
            prefix="cr-relativistic-writer-adversarial-") as tmp:
        yield Path(tmp)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--input", default=DEFAULT_INPUT, type=Path)
    parser.add_argument("--work-dir", type=Path)
    args = parser.parse_args()

    _require(args.binary.is_file(), f"missing Athena binary: {args.binary}")
    _require(args.input.is_file(), f"missing input file: {args.input}")
    with _runtime_root(args.work_dir) as root:
        _pmom_scaled_norm_success(args.binary, root, args.input)
        _pass("pmom scaled-norm success")
        _pmom_preaccumulation_rejection(args.binary, root, args.input)
        _pass("pmom pre-accumulation pre-artifact rejection")
        _df_scaled_norm_success(args.binary, root, args.input)
        _pass("df scaled-norm success")
        _pspec_mu_scaled_speed_success(args.binary, root, args.input)
        _pass("pspec mu scaled-speed success")
        _pspec_proxy_stability_success(args.binary, root, args.input)
        _pass("pspec magnetic-moment proxy stability success")
        _pspec_proxy_tiny_b_success(args.binary, root, args.input)
        _pass("pspec tiny-B proxy stability success")
        _pspec2_scaled_speed_success(args.binary, root, args.input)
        _pass("pspec2 scaled-speed success")
        _pspec2_vperp_stability_success(args.binary, root, args.input)
        _pass("pspec2 vperp stability success")
        _psamp_vperp_stability_success(args.binary, root, args.input)
        _pass("psamp vperp stability success")
        _psamp_proxy_tiny_b_success(args.binary, root, args.input)
        _pass("psamp tiny-B proxy stability success")
        _psamp_row_atomic_rejection(args.binary, root, args.input)
        _pass("psamp row-atomic pre-artifact rejection")
        _psamp_ratio_scaled_success(args.binary, root, args.input)
        _pass("psamp applicability-ratio scaled success")
        _zero_field_rejections(args.binary, root, args.input)
        _ppd_preopen_rejection(args.binary, root, args.input)
        _pspec_log_floor_metadata_success(args.binary, root, args.input)
        _pass(
            "ppd float32 species predicate compile-bound",
            "static_assert exact=16777216 lossy_rejected=16777217 exact=2147483520",
        )

    print("PASS: retained relativistic writer adversarial probes (20/20)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
