#!/usr/bin/env python3
"""Reject deliberate Phase-7 relativistic diagnostic metadata corruptions."""

from __future__ import annotations

import json
import struct
import sys
import tempfile
from pathlib import Path
from typing import Callable

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import cr_relativistic_diagnostics_runtime_inspect as diagnostics  # noqa: E402
import cr_relativistic_all_formats_inspect as all_formats  # noqa: E402
import cr_tracer_inspect as inspector  # noqa: E402


REL_COMMON = (
    "schema=akcr_particle_output_v1 mode=relativistic_hc units=code_model "
    "c_model=1 alpha_s=1"
)
PPD_METADATA = (
    f"{REL_COMMON} columns=species_x_y_z position_units=code_length payload=float32"
)
PMOM_METADATA = (
    f"{REL_COMMON} quantity_basis=physical_velocity_shadow mu_units=dimensionless "
    "displacement_units=code_length "
    "displacement_second_moment_units=code_length_squared "
    "velocity_units=code_velocity velocity_second_moment_units=code_velocity_squared "
    "dparallel_definition=accumulated_midpoint_sum_dx_dot_Bhat"
)
PSAMP_METADATA = (
    f"{REL_COMMON} selected_species=-1 sample_stride=1 sample_offset=0 "
    "sample_count=1 field_count=1 position_units=code_length "
    "velocity_units=code_velocity w_units=code_velocity "
    "cE_units=code_velocity_times_code_B gamma_units=dimensionless "
    "kinetic_energy_model_units=code_velocity_squared "
    "work_units=code_velocity_squared alpha_s_units=inverse_code_time_per_code_B "
    "r_larmor_over_dx_min_units=dimensionless"
)


def _require_rejection(label: str, action: Callable[[], object],
                       expected: str) -> None:
    try:
        action()
    except (RuntimeError, ValueError) as error:
        if expected not in str(error):
            raise RuntimeError(
                f"{label}: rejection {error!r} lacks expected substring {expected!r}"
            ) from error
        print(f"rejected {label}: {error}")
        return
    raise RuntimeError(f"{label}: mutation was accepted")


def _write_sample(path: Path, metadata: str, columns: str, row: str,
                  nranks: int = 1, nranks_text: str | None = None,
                  nranks_key: str = "nranks") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header_nranks = str(nranks) if nranks_text is None else nranks_text
    path.write_text(
        f"# AthenaK particle sample at time= 0 {nranks_key}= {header_nranks} cycle= 1\n"
        f"# metadata: {metadata}\n"
        f"# columns: {columns}\n"
        f"{row}\n"
    )


def _write_spectrum(path: Path, metadata: str, values: tuple[int, ...],
                    layout: str | None = "species-major int32 histogram") -> None:
    header = (
        "# AthenaK particle spectrum at time= 0 nranks= 1 cycle= 1\n"
        f"# metadata: {metadata}\n"
        + (f"# layout: {layout}\n" if layout is not None else "")
    ).encode("ascii")
    path.write_bytes(header + struct.pack("<" + str(len(values)) + "i", *values))


def _write_joint_spectrum(path: Path, metadata: str,
                          values: tuple[int, ...],
                          layout: str | None = (
                              "species-major int32 histogram, bin1 then bin2"
                          )) -> None:
    header = (
        "# AthenaK particle joint spectrum at time= 0 nranks= 1 cycle= 1\n"
        f"# metadata: {metadata}\n"
        + (f"# layout: {layout}\n" if layout is not None else "")
    ).encode("ascii")
    path.write_bytes(header + struct.pack("<" + str(len(values)) + "i", *values))


def _write_histogram(path: Path, marker: str, metadata: str,
                     values: tuple[int, ...],
                     layout: str | None = "default") -> None:
    if layout == "default":
        layout = {
            "# AthenaK particle distribution function":
                "species-major int32 histogram; mu=dot(v,B)/(|v||B|)",
            "# AthenaK particle displacement histogram":
                "species-major int32 histograms ordered dx, dy, dz",
        }.get(marker, "species-major int32 histogram")
    header = (
        f"{marker} at time= 0 nranks= 1 cycle= 1\n"
        f"# metadata: {metadata}\n"
        + (f"# layout: {layout}\n" if layout is not None else "")
    ).encode("ascii")
    path.write_bytes(header + struct.pack("<" + str(len(values)) + "i", *values))


def _write_ppd(path: Path, metadata: str, values: tuple[float, ...],
               particles: int = 1, particles_text: str | None = None,
               particles_key: str = "particles") -> None:
    header_particles = str(particles) if particles_text is None else particles_text
    header = (
        "# AthenaK particle position data at time= 0 nranks= 1 "
        f"{particles_key}= {header_particles} cycle= 1\n"
        f"# metadata: {metadata}\n"
    ).encode("ascii")
    path.write_bytes(header + struct.pack("<" + str(len(values)) + "f", *values))


def _write_moments(path: Path, rows: tuple[str, ...],
                   include_header: bool = True, metadata: str | None = None) -> None:
    prefix = (
        "# AthenaK particle moments at time= 0 nranks= 1 cycle= 1\n"
        if include_header else ""
    )
    if metadata is not None:
        prefix += f"# metadata: {metadata}\n"
    path.write_text(prefix + "\n".join(rows) + "\n")


def _modern_moment_row(species: int, count: int = 1) -> str:
    return " ".join(str(value) for value in (species, count, *([0] * 24)))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="cr-rel-phase7-negative-") as temporary:
        root = Path(temporary)
        sample = root / "sample.psamp"
        modern_sample = root / "psamp/rank_00000000/sample.psamp"
        spectrum = root / "sample.pspec"
        joint_spectrum = root / "sample.pspec2"
        ppd = root / "sample.ppd"
        df = root / "sample.df"
        dxh = root / "sample.dxh"
        drh = root / "sample.drh"
        dparh = root / "sample.dparh"
        pmom = root / "sample.pmom"

        criteria = json.loads(diagnostics.DEFAULT_CRITERIA.read_text())
        criteria["criteria"][0]["limit"] = 1.0
        tampered = root / "tampered-criteria.json"
        tampered.write_text(json.dumps(criteria, indent=2) + "\n")
        _require_rejection(
            "tampered criteria digest",
            lambda: diagnostics._validate_criteria(tampered, {}),
            "criteria SHA-256",
        )

        provenance_file = root / "provenance.bin"
        provenance_file.write_bytes(b"phase7 provenance")
        _require_rejection(
            "binary provenance digest mismatch",
            lambda: all_formats._require_matching_sha256(
                {"binary_sha256": "0" * 64}, "binary_sha256", provenance_file),
            "metrics binary_sha256",
        )
        _require_rejection(
            "input provenance digest mismatch",
            lambda: all_formats._require_matching_sha256(
                {"input_sha256": "f" * 64}, "input_sha256", provenance_file),
            "metrics input_sha256",
        )

        _write_sample(
            sample, "sample_count=1 field_count=2",
            "rank pgid tag species x x1", "0 0 1 0 0.1 0.1")
        _require_rejection(
            "canonical sample alias collision",
            lambda: inspector.read_psamp_file(sample),
            "duplicate canonical sample columns",
        )

        _write_sample(
            sample, "sample_count=1 field_count=2",
            "rank pgid tag species x", "0 0 1 0 0.1")
        _require_rejection(
            "sample field-count mismatch",
            lambda: inspector.read_psamp_file(sample),
            "field_count=2",
        )

        _write_sample(
            sample, "sample_count=2 field_count=1",
            "rank pgid tag species x", "0 0 1 0 0.1")
        _require_rejection(
            "sample row-count mismatch",
            lambda: inspector.read_psamp_file(sample),
            "sample_count=2",
        )

        _write_sample(
            sample, "sample_count=1 field_count=1",
            "rank pgid tag species mystery", "0 0 1 0 0.1")
        _require_rejection(
            "unknown sample column",
            lambda: inspector.read_psamp_file(sample),
            "unknown sample column",
        )

        _write_sample(
            sample, "mode=relativistic_hc sample_count=1 field_count=1",
            "rank pgid tag species p", "0 0 1 0 0.1")
        _require_rejection(
            "relativistic sample legacy alias",
            lambda: inspector.read_psamp_file(sample),
            "relativistic_hc sample rejects ambiguous",
        )

        sample.write_text(
            "# AthenaK particle sample at time= 0 nranks= 1 cycle= 1\n"
            "# metadata: sample_count=1 field_count=1\n"
            "# columns: rank pgid tag species p\n"
            "# metadata: mode=relativistic_hc\n"
            "0 0 1 0 0.1\n"
        )
        _require_rejection(
            "sample metadata after columns",
            lambda: inspector.read_psamp_file(sample),
            "metadata after columns",
        )

        _write_sample(
            sample, "sample_count=1 field_count=1",
            "rank pgid tag species x", "0 0 1 0 nan")
        _require_rejection(
            "non-finite sample value",
            lambda: inspector.read_psamp_file(sample),
            "non-finite x",
        )

        _write_sample(
            sample, "sample_count=1 sample_count=1 field_count=1",
            "rank pgid tag species x", "0 0 1 0 0.1")
        _require_rejection(
            "duplicate sample metadata key",
            lambda: inspector.read_psamp_file(sample),
            "duplicate metadata key",
        )

        _write_spectrum(spectrum, "quantity=wmag nbin=2 reduce=1", (1, 0))
        record = inspector.read_pspec_record(spectrum, 1, 2)
        if record["metadata"].get("quantity") != "wmag":
            raise RuntimeError("spectrum metadata retention failed")
        print("accepted spectrum metadata retention")

        _write_spectrum(spectrum, "quantity=wmag nbin=2 reduce=1", (1, 0, 0))
        _require_rejection(
            "spectrum nbin mismatch",
            lambda: inspector.read_pspec_record(spectrum, 1, 3),
            "metadata nbin=2",
        )

        _write_spectrum(
            spectrum, "mode=relativistic_hc quantity=p nbin=2 reduce=1", (1, 0))
        _require_rejection(
            "relativistic spectrum legacy alias",
            lambda: inspector.read_pspec_record(spectrum, 1, 2),
            "relativistic_hc spectrum rejects ambiguous",
        )

        _write_spectrum(
            spectrum,
            "mode=relativistic_hc quantity=wmag nspecies=2 nbin=2 reduce=1",
            (1, 0),
        )
        _require_rejection(
            "relativistic spectrum species metadata mismatch",
            lambda: inspector.read_pspec_record(spectrum, 1, 2),
            "metadata nspecies=2",
        )

        _write_joint_spectrum(
            joint_spectrum,
            "mode=relativistic_hc quantity=mu_p nbin1=1 nbin2=2 reduce=1",
            (1, 0),
        )
        _require_rejection(
            "relativistic joint-spectrum legacy alias",
            lambda: inspector.read_pspec2_record(joint_spectrum, 1, 1, 2),
            "relativistic_hc joint spectrum rejects ambiguous",
        )

        _write_joint_spectrum(
            joint_spectrum,
            "mode=relativistic_hc quantity=mu_wmag nspecies=2 "
            "nbin1=1 nbin2=2 reduce=1",
            (1, 0),
        )
        _require_rejection(
            "relativistic joint-spectrum species metadata mismatch",
            lambda: inspector.read_pspec2_record(joint_spectrum, 1, 1, 2),
            "metadata nspecies=2",
        )

        _write_ppd(ppd, PPD_METADATA, (0.0, 0.1, 0.2, 0.3))
        ppd_record = inspector.read_ppd_file(ppd)
        if ppd_record["metadata"].get("mode") != "relativistic_hc":
            raise RuntimeError("ppd metadata retention failed")
        print("accepted ppd metadata retention")

        _write_ppd(ppd, "mode=relativistic_hc mode=legacy", (0.0, 0.1, 0.2, 0.3))
        _require_rejection(
            "duplicate ppd metadata key",
            lambda: inspector.read_ppd_file(ppd),
            "duplicate metadata key",
        )

        _write_ppd(
            ppd, PPD_METADATA, (0.0, 0.1, 0.2, 0.3), particles=2)
        _require_rejection(
            "ppd header payload count mismatch",
            lambda: inspector.read_ppd_file(ppd),
            "ppd header particles=2, payload count=1",
        )

        _write_ppd(ppd, PPD_METADATA, (0.25, 0.1, 0.2, 0.3))
        _require_rejection(
            "fractional ppd species",
            lambda: inspector.read_ppd_file(ppd),
            "non-integral, negative, or out-of-range ppd species",
        )

        for label, particles_text, expected in (
            ("duplicate ppd particle count", "1 particles= 2",
             "header requires exactly one particles integer token"),
            ("fractional ppd particle count", "1.5",
             "header has invalid particles='1.5'"),
        ):
            _write_ppd(
                ppd, PPD_METADATA, (0.0, 0.1, 0.2, 0.3),
                particles_text=particles_text)
            _require_rejection(
                label,
                lambda: inspector.read_ppd_file(ppd),
                expected,
            )

        _write_ppd(
            ppd, PPD_METADATA, (0.0, 0.1, 0.2, 0.3),
            particles_key="x-particles")
        _require_rejection(
            "prefixed ppd particle-count label",
            lambda: inspector.read_ppd_file(ppd),
            "header requires exactly one particles integer token",
        )

        for label, particles_text, expected in (
            ("missing ppd particle count", "1", "exactly one particles"),
            ("out-of-range ppd particle count", "2147483648",
             "header has out-of-range particles='2147483648'"),
        ):
            _write_ppd(
                ppd, PPD_METADATA, (0.0, 0.1, 0.2, 0.3),
                particles_text=particles_text,
                particles_key=(
                    "other" if label == "missing ppd particle count" else "particles"),
            )
            _require_rejection(
                label,
                lambda: inspector.read_ppd_file(ppd),
                expected,
            )

        _write_ppd(ppd, PPD_METADATA, (2147483520.0, 0.1, 0.2, 0.3))
        ppd_boundary = inspector.read_ppd_file(ppd)
        if ppd_boundary["particles"][0]["species"] != 2147483520:
            raise RuntimeError("exact float32 ppd species boundary changed")
        print("accepted exact float32 ppd species boundary")

        _write_ppd(ppd, PPD_METADATA, (4294967296.0, 0.1, 0.2, 0.3))
        _require_rejection(
            "out-of-range ppd species",
            lambda: inspector.read_ppd_file(ppd),
            "non-integral, negative, or out-of-range ppd species",
        )

        _write_ppd(
            ppd,
            PPD_METADATA.replace(
                "schema=akcr_particle_output_v1",
                "schema=akcr_particle_output_v1_typo",
            ).replace("mode=relativistic_hc", "mode=relativistic_hc_typo"),
            (0.0, 0.1, 0.2, 0.3),
        )
        _require_rejection(
            "malformed reserved ppd markers",
            lambda: inspector.read_ppd_file(ppd),
            "metadata schema='akcr_particle_output_v1_typo', expected "
            "'akcr_particle_output_v1'",
        )

        valid_moment_row = "0 1 0 0 0 0 0 0 0 0 0 0 0 0"
        _write_moments(pmom, (valid_moment_row,), include_header=False)
        _require_rejection(
            "pmom missing header",
            lambda: inspector.read_pmom_record(pmom, 1),
            "missing particle-moment header",
        )

        _write_moments(pmom, ("0 1.25 0 0 0 0 0 0 0 0 0 0 0 0",))
        _require_rejection(
            "fractional pmom count",
            lambda: inspector.read_pmom_record(pmom, 1),
            "non-integral or negative count",
        )

        _write_moments(pmom, ("0.25 1 0 0 0 0 0 0 0 0 0 0 0 0",))
        _require_rejection(
            "fractional pmom species",
            lambda: inspector.read_pmom_record(pmom, 1),
            "non-integral or negative species",
        )

        _write_moments(pmom, (valid_moment_row, valid_moment_row))
        _require_rejection(
            "pmom latest block row count mismatch",
            lambda: inspector.read_pmom_record(pmom, 1),
            "latest block contains 2 particle-moment rows, expected 1",
        )

        for label, rows in (
            ("duplicate", (_modern_moment_row(0), _modern_moment_row(0))),
            ("reordered", (_modern_moment_row(1), _modern_moment_row(0))),
            ("out-of-range", (_modern_moment_row(0), _modern_moment_row(2))),
        ):
            _write_moments(pmom, rows, metadata=PMOM_METADATA)
            _require_rejection(
                f"{label} relativistic pmom species identity",
                lambda: inspector.read_pmom_record(pmom, 2),
                "relativistic_hc pmom row",
            )

        _write_histogram(
            df, "# AthenaK particle distribution function",
            "quantity=mu nspecies=1 nbin=2 vmin=-1 vmax=1 reduce=1", (35, 0))
        df_record = inspector.read_df_record(df, 1, 2)
        if df_record["histogram"] != [[35, 0]]:
            raise RuntimeError("leading-count-35 histogram parsing regression")
        print("accepted histogram leading-count-35 binary payload")

        for first_bin in (32, 10, 13):
            _write_histogram(
                df, "# AthenaK particle distribution function",
                "quantity=mu nspecies=1 nbin=2 vmin=-1 vmax=1 reduce=1",
                (first_bin, 0),
            )
            df_record = inspector.read_df_record(df, 1, 2)
            if df_record["histogram"] != [[first_bin, 0]]:
                raise RuntimeError(
                    f"leading-count-{first_bin} histogram parsing regression")
            print(f"accepted histogram leading-count-{first_bin} binary payload")

        _write_histogram(
            df, "# AthenaK particle distribution function",
            "quantity=mu nspecies=2 nbin=2 vmin=-1 vmax=1 reduce=1", (35, 0))
        _require_rejection(
            "histogram species metadata mismatch",
            lambda: inspector.read_df_record(df, 1, 2),
            "metadata nspecies=2",
        )

        negative_histograms = (
            (
                "df",
                lambda: _write_histogram(
                    df, "# AthenaK particle distribution function",
                    "quantity=mu nspecies=1 nbin=2", (-1, 2)),
                lambda: inspector.read_df_record(df, 1, 2),
            ),
            (
                "dxh",
                lambda: _write_histogram(
                    dxh, "# AthenaK particle displacement histogram",
                    "quantity=dx_dy_dz nspecies=1 nbin=2", (-1, 2, 0, 0, 0, 0)),
                lambda: inspector.read_dxh_record(dxh, 1, 2),
            ),
            (
                "drh",
                lambda: _write_histogram(
                    drh, "# AthenaK particle scalar displacement histogram",
                    "quantity=displacement_norm nspecies=1 nbin=2", (-1, 2)),
                lambda: inspector.read_scalar_histogram_record(
                    drh, b"# AthenaK particle scalar displacement histogram", 1, 2),
            ),
            (
                "dparh",
                lambda: _write_histogram(
                    dparh, "# AthenaK particle parallel displacement histogram",
                    "quantity=dparallel nspecies=1 nbin=2", (-1, 2)),
                lambda: inspector.read_scalar_histogram_record(
                    dparh, b"# AthenaK particle parallel displacement histogram",
                    1, 2),
            ),
            (
                "pspec",
                lambda: _write_spectrum(
                    spectrum, "quantity=wmag nspecies=1 nbin=2", (-1, 2)),
                lambda: inspector.read_pspec_record(spectrum, 1, 2),
            ),
            (
                "pspec2",
                lambda: _write_joint_spectrum(
                    joint_spectrum,
                    "quantity=mu_wmag nspecies=1 nbin1=1 nbin2=2",
                    (-1, 2),
                ),
                lambda: inspector.read_pspec2_record(joint_spectrum, 1, 1, 2),
            ),
        )
        for name, write_mutation, read_mutation in negative_histograms:
            write_mutation()
            _require_rejection(
                f"negative {name} histogram bin",
                read_mutation,
                "negative histogram bin count",
            )

        invalid_layout_histograms = (
            (
                "df",
                lambda: _write_histogram(
                    df, "# AthenaK particle distribution function",
                    "quantity=mu nspecies=1 nbin=2", (1, 0), layout="wrong"),
                lambda: inspector.read_df_record(df, 1, 2),
            ),
            (
                "dxh",
                lambda: _write_histogram(
                    dxh, "# AthenaK particle displacement histogram",
                    "quantity=dx_dy_dz nspecies=1 nbin=2",
                    (1, 0, 1, 0, 1, 0), layout="wrong"),
                lambda: inspector.read_dxh_record(dxh, 1, 2),
            ),
            (
                "drh",
                lambda: _write_histogram(
                    drh, "# AthenaK particle scalar displacement histogram",
                    "quantity=displacement_norm nspecies=1 nbin=2",
                    (1, 0), layout="wrong"),
                lambda: inspector.read_scalar_histogram_record(
                    drh, b"# AthenaK particle scalar displacement histogram", 1, 2),
            ),
            (
                "dparh",
                lambda: _write_histogram(
                    dparh, "# AthenaK particle parallel displacement histogram",
                    "quantity=dparallel nspecies=1 nbin=2",
                    (1, 0), layout="wrong"),
                lambda: inspector.read_scalar_histogram_record(
                    dparh, b"# AthenaK particle parallel displacement histogram",
                    1, 2),
            ),
            (
                "pspec",
                lambda: _write_spectrum(
                    spectrum, "quantity=wmag nspecies=1 nbin=2",
                    (1, 0), layout="wrong"),
                lambda: inspector.read_pspec_record(spectrum, 1, 2),
            ),
            (
                "pspec2",
                lambda: _write_joint_spectrum(
                    joint_spectrum,
                    "quantity=mu_wmag nspecies=1 nbin1=1 nbin2=2",
                    (1, 0), layout="wrong"),
                lambda: inspector.read_pspec2_record(joint_spectrum, 1, 1, 2),
            ),
        )
        for name, write_mutation, read_mutation in invalid_layout_histograms:
            write_mutation()
            _require_rejection(
                f"invalid {name} histogram layout",
                read_mutation,
                "histogram layout 'wrong'",
            )

        _write_histogram(
            df, "# AthenaK particle distribution function",
            "mode=relativistic_hc quantity=mu nspecies=1 nbin=2",
            (1, 0), layout=None)
        _require_rejection(
            "missing relativistic histogram layout",
            lambda: inspector.read_df_record(df, 1, 2),
            "relativistic_hc histogram is missing layout metadata",
        )

        _write_ppd(ppd, "mode=legacy junk", (0.0, 0.1, 0.2, 0.3))
        _require_rejection(
            "malformed ppd metadata token",
            lambda: inspector.read_ppd_file(ppd),
            "Malformed metadata token 'junk'",
        )

        _write_ppd(ppd, "mode=relativistic_hc", (0.0, 0.1, 0.2, 0.3))
        _require_rejection(
            "underspecified relativistic ppd metadata",
            lambda: inspector.read_ppd_file(ppd),
            "relativistic_hc ppd metadata is missing",
        )

        underspecified_histograms = (
            (
                "df",
                lambda: _write_histogram(
                    df, "# AthenaK particle distribution function",
                    "mode=relativistic_hc", (1, 0)),
                lambda: inspector.read_df_record(df, 1, 2),
            ),
            (
                "dxh",
                lambda: _write_histogram(
                    dxh, "# AthenaK particle displacement histogram",
                    "mode=relativistic_hc", (1, 0, 1, 0, 1, 0)),
                lambda: inspector.read_dxh_record(dxh, 1, 2),
            ),
            (
                "drh",
                lambda: _write_histogram(
                    drh, "# AthenaK particle scalar displacement histogram",
                    "mode=relativistic_hc", (1, 0)),
                lambda: inspector.read_scalar_histogram_record(
                    drh, b"# AthenaK particle scalar displacement histogram", 1, 2),
            ),
            (
                "dparh",
                lambda: _write_histogram(
                    dparh, "# AthenaK particle parallel displacement histogram",
                    "mode=relativistic_hc", (1, 0)),
                lambda: inspector.read_scalar_histogram_record(
                    dparh, b"# AthenaK particle parallel displacement histogram",
                    1, 2),
            ),
            (
                "pspec",
                lambda: _write_spectrum(
                    spectrum, "mode=relativistic_hc quantity=wmag", (1, 0)),
                lambda: inspector.read_pspec_record(spectrum, 1, 2),
            ),
            (
                "pspec2",
                lambda: _write_joint_spectrum(
                    joint_spectrum, "mode=relativistic_hc quantity=mu_wmag", (1, 0)),
                lambda: inspector.read_pspec2_record(joint_spectrum, 1, 1, 2),
            ),
        )
        for name, write_mutation, read_mutation in underspecified_histograms:
            write_mutation()
            _require_rejection(
                f"underspecified relativistic {name} metadata",
                read_mutation,
                f"relativistic_hc {name} metadata is missing",
            )

        _write_moments(pmom, (valid_moment_row,), metadata=PMOM_METADATA)
        _require_rejection(
            "relativistic pmom legacy-width row",
            lambda: inspector.read_pmom_record(pmom, 1),
            "relativistic_hc pmom requires 26-column rows",
        )

        _write_sample(
            sample, "mode=relativistic_hc sample_count=1 field_count=1",
            "rank pgid tag species x", "0 0 1 0 0.1")
        _require_rejection(
            "underspecified relativistic psamp metadata",
            lambda: inspector.read_psamp_file(sample),
            "relativistic_hc psamp metadata is missing",
        )

        for label, row, expected in (
            ("negative sample rank", "-1 0 1 0 0.1", "out-of-range rank"),
            ("negative sample pgid", "0 -1 1 0 0.1", "out-of-range pgid"),
            ("negative sample tag", "0 0 -1 0 0.1", "out-of-range tag"),
            ("negative sample species", "0 0 1 -1 0.1", "out-of-range species"),
            (
                "sample shard rank mismatch",
                "1 0 1 0 0.1",
                "does not match shard rank 0",
            ),
        ):
            _write_sample(
                modern_sample, PSAMP_METADATA, "rank pgid tag species x", row)
            _require_rejection(
                label,
                lambda: inspector.read_psamp_file(modern_sample),
                expected,
            )

        _write_sample(
            sample, PSAMP_METADATA, "rank pgid tag species x", "0 0 1 0 0.1")
        _require_rejection(
            "unsharded relativistic sample",
            lambda: inspector.read_psamp_file(sample),
            "requires canonical psamp/rank_######## shard path",
        )

        substring_sample = root / "foo_rank_00000007_bar/sample.psamp"
        _write_sample(
            substring_sample, PSAMP_METADATA,
            "rank pgid tag species x", "7 0 1 0 0.1", nranks=8)
        _require_rejection(
            "substring relativistic sample shard",
            lambda: inspector.read_psamp_file(substring_sample),
            "requires canonical psamp/rank_######## shard path",
        )

        rank7_sample = root / "psamp/rank_00000007/sample.psamp"
        _write_sample(
            rank7_sample, PSAMP_METADATA,
            "rank pgid tag species x", "7 0 1 0 0.1")
        _require_rejection(
            "relativistic sample shard rank outside header",
            lambda: inspector.read_psamp_file(rank7_sample),
            "sample shard rank 7 is outside nranks=1",
        )

        for label, nranks_text, expected in (
            ("duplicate sample rank count", "1 nranks= 2",
             "header requires exactly one nranks integer token"),
            ("fractional sample rank count", "1.5",
             "header has invalid nranks='1.5'"),
        ):
            _write_sample(
                modern_sample, PSAMP_METADATA, "rank pgid tag species x",
                "0 0 1 0 0.1", nranks_text=nranks_text)
            _require_rejection(
                label,
                lambda: inspector.read_psamp_file(modern_sample),
                expected,
            )

        _write_sample(
            modern_sample, PSAMP_METADATA, "rank pgid tag species x",
            "0 0 1 0 0.1", nranks_key="x-nranks")
        _require_rejection(
            "prefixed sample rank-count label",
            lambda: inspector.read_psamp_file(modern_sample),
            "sample header is missing nranks",
        )

        for label, nranks_text, nranks_key, expected in (
            ("missing sample rank count", "1", "other",
             "sample header is missing nranks"),
            ("zero sample rank count", "0", "nranks",
             "header has out-of-range nranks='0'"),
            ("out-of-range sample rank count", "2147483648", "nranks",
             "header has out-of-range nranks='2147483648'"),
        ):
            _write_sample(
                modern_sample, PSAMP_METADATA, "rank pgid tag species x",
                "0 0 1 0 0.1", nranks_text=nranks_text, nranks_key=nranks_key)
            _require_rejection(
                label,
                lambda: inspector.read_psamp_file(modern_sample),
                expected,
            )

        for label, metadata, row, expected in (
            (
                "zero sample stride",
                PSAMP_METADATA.replace("sample_stride=1", "sample_stride=0"),
                "0 0 1 0 0.1",
                "out-of-range sample_stride",
            ),
            (
                "sample offset equals stride",
                PSAMP_METADATA.replace("sample_offset=0", "sample_offset=1"),
                "0 0 1 0 0.1",
                "0 <= sample_offset < sample_stride",
            ),
            (
                "sample selected-species mismatch",
                PSAMP_METADATA.replace("selected_species=-1", "selected_species=1"),
                "0 0 1 0 0.1",
                "does not satisfy selection metadata",
            ),
            (
                "sample splitmix64 mismatch",
                PSAMP_METADATA.replace(
                    "sample_stride=1 sample_offset=0",
                    f"sample_stride=2 sample_offset="
                    f"{1 - inspector.sample_hash(0, 1) % 2}",
                ),
                "0 0 1 0 0.1",
                "does not satisfy selection metadata",
            ),
            (
                "sample rank upper bound",
                PSAMP_METADATA,
                "2147483648 0 1 0 0.1",
                "out-of-range rank",
            ),
            (
                "sample pgid upper bound",
                PSAMP_METADATA,
                "0 2147483648 1 0 0.1",
                "out-of-range pgid",
            ),
            (
                "sample tag upper bound",
                PSAMP_METADATA,
                "0 0 2147483648 0 0.1",
                "out-of-range tag",
            ),
            (
                "sample species upper bound",
                PSAMP_METADATA,
                "0 0 1 2147483648 0.1",
                "out-of-range species",
            ),
        ):
            _write_sample(
                modern_sample, metadata, "rank pgid tag species x", row)
            _require_rejection(
                label,
                lambda: inspector.read_psamp_file(modern_sample),
                expected,
            )

        _write_spectrum(
            spectrum,
            f"{REL_COMMON} quantity=wmag quantity_units=dimensionless "
            "histogram_units=particle_count nspecies=1 nbin=2 "
            "vmin=0 vmax=1 reduce=1",
            (1, 0),
        )
        _require_rejection(
            "relativistic spectrum contradictory quantity units",
            lambda: inspector.read_pspec_record(spectrum, 1, 2),
            "quantity_units='dimensionless', expected 'code_velocity'",
        )

        typo_rel_common = REL_COMMON.replace(
            "schema=akcr_particle_output_v1",
            "schema=akcr_particle_output_v1_typo",
        ).replace("mode=relativistic_hc", "mode=relativistic_hc_typo")
        _write_spectrum(
            spectrum,
            f"{typo_rel_common} quantity=wmag quantity_units=code_velocity "
            "histogram_units=particle_count nspecies=1 nbin=2 "
            "vmin=0 vmax=1 reduce=1",
            (1, 0),
        )
        _require_rejection(
            "malformed reserved spectrum markers",
            lambda: inspector.read_pspec_record(spectrum, 1, 2),
            "metadata schema='akcr_particle_output_v1_typo', expected "
            "'akcr_particle_output_v1'",
        )

        _write_joint_spectrum(
            joint_spectrum,
            f"{REL_COMMON} quantity=vpar_vperp axis1_units=code_velocity "
            "axis2_units=bananas histogram_units=particle_count nspecies=1 "
            "nbin1=1 nbin2=2 vmin1=-1 vmax1=1 vmin2=0 vmax2=1 reduce=1",
            (1, 0),
        )
        _require_rejection(
            "relativistic joint-spectrum contradictory axis units",
            lambda: inspector.read_pspec2_record(joint_spectrum, 1, 1, 2),
            "axis2_units='bananas', expected 'code_velocity'",
        )

        _write_spectrum(
            spectrum,
            f"{REL_COMMON} quantity=log10_kinetic_energy_model "
            "quantity_units=log10_code_velocity_squared "
            "histogram_units=particle_count nspecies=1 nbin=2 "
            "vmin=-4 vmax=1 reduce=1",
            (1, 0),
        )
        _require_rejection(
            "relativistic log spectrum missing floor",
            lambda: inspector.read_pspec_record(spectrum, 1, 2),
            "relativistic_hc pspec metadata is missing "
            "['log10_floor', 'log10_floor_units']",
        )

        for label, floor_metadata, expected in (
            (
                "relativistic log spectrum missing floor units",
                "log10_floor=1e-300",
                "relativistic_hc pspec metadata is missing ['log10_floor_units']",
            ),
            (
                "relativistic log spectrum contradictory floor units",
                "log10_floor=1e-300 log10_floor_units=dimensionless",
                "log10_floor_units='dimensionless', expected "
                "'code_velocity_squared'",
            ),
            (
                "relativistic log spectrum zero floor",
                "log10_floor=0 log10_floor_units=code_velocity_squared",
                "metadata requires log10_floor > 0",
            ),
            (
                "relativistic log spectrum negative floor",
                "log10_floor=-1 log10_floor_units=code_velocity_squared",
                "metadata requires log10_floor > 0",
            ),
            (
                "relativistic log spectrum non-finite floor",
                "log10_floor=nan log10_floor_units=code_velocity_squared",
                "metadata has non-finite log10_floor='nan'",
            ),
        ):
            _write_spectrum(
                spectrum,
                f"{REL_COMMON} quantity=log10_kinetic_energy_model "
                "quantity_units=log10_code_velocity_squared "
                "histogram_units=particle_count nspecies=1 nbin=2 "
                f"vmin=-400 vmax=1 reduce=1 {floor_metadata}",
                (1, 0),
            )
            _require_rejection(
                label,
                lambda: inspector.read_pspec_record(spectrum, 1, 2),
                expected,
            )

        legacy_no_layout_histograms = (
            (
                "df",
                lambda: _write_histogram(
                    df, "# AthenaK particle distribution function",
                    "quantity=mu nspecies=1 nbin=2", (1, 0), layout=None),
                lambda: inspector.read_df_record(df, 1, 2),
            ),
            (
                "dxh",
                lambda: _write_histogram(
                    dxh, "# AthenaK particle displacement histogram",
                    "quantity=dx_dy_dz nspecies=1 nbin=2",
                    (1, 0, 1, 0, 1, 0), layout=None),
                lambda: inspector.read_dxh_record(dxh, 1, 2),
            ),
            (
                "drh",
                lambda: _write_histogram(
                    drh, "# AthenaK particle scalar displacement histogram",
                    "quantity=displacement_norm nspecies=1 nbin=2",
                    (1, 0), layout=None),
                lambda: inspector.read_scalar_histogram_record(
                    drh, b"# AthenaK particle scalar displacement histogram", 1, 2),
            ),
            (
                "dparh",
                lambda: _write_histogram(
                    dparh, "# AthenaK particle parallel displacement histogram",
                    "quantity=dparallel nspecies=1 nbin=2", (1, 0), layout=None),
                lambda: inspector.read_scalar_histogram_record(
                    dparh, b"# AthenaK particle parallel displacement histogram",
                    1, 2),
            ),
            (
                "pspec",
                lambda: _write_spectrum(
                    spectrum, "quantity=wmag nspecies=1 nbin=2",
                    (1, 0), layout=None),
                lambda: inspector.read_pspec_record(spectrum, 1, 2),
            ),
            (
                "pspec2",
                lambda: _write_joint_spectrum(
                    joint_spectrum, "quantity=mu_wmag nspecies=1 nbin1=1 nbin2=2",
                    (1, 0), layout=None),
                lambda: inspector.read_pspec2_record(joint_spectrum, 1, 1, 2),
            ),
        )
        for name, write_legacy, read_legacy in legacy_no_layout_histograms:
            write_legacy()
            read_legacy()
            print(f"accepted legacy {name} histogram without layout metadata")

    print("CR relativistic Phase-7 diagnostic mutation controls passed: 88 rejected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
