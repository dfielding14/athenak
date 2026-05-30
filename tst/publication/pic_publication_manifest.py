"""Manifest for publication-oriented PIC/MHD-PIC runs and artifacts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
from typing import Tuple


@dataclass(frozen=True)
class PublicationCase:
    case_id: str
    group: str
    runner: str  # "module" or "deck"
    module: str = ""
    input_deck: str = ""
    output_globs: Tuple[str, ...] = ()
    ranks_local: int = 1
    ranks_hpc: int = 8
    athena_args: Tuple[str, ...] = ()
    profile: str = "local"  # "local", "hpc", or "both"
    note: str = ""


PUBLICATION_CASES: Tuple[PublicationCase, ...] = (
    # ---------------------------------------------------------------------
    # Proxy-regression anchors (unchanged physics and thresholds)
    # ---------------------------------------------------------------------
    PublicationCase(
        case_id="entity_deposit_mink",
        group="entity_core_proxy",
        runner="module",
        module="scripts.particles.pic_entity_deposit_mink",
        output_globs=("tst/build/src/bin/pic_entity_dep_mink_*.*.bin",),
        note="Entity parity: Minkowski deposit invariance.",
    ),
    PublicationCase(
        case_id="entity_deposit_reflect",
        group="entity_core_proxy",
        runner="module",
        module="scripts.particles.pic_entity_deposit_reflect",
        output_globs=("tst/build/src/bin/pic_entity_dep_reflect_*.*.bin",),
        note="Entity parity: reflecting-boundary deposit invariance.",
    ),
    PublicationCase(
        case_id="em_vacuum_wave",
        group="entity_core_proxy",
        runner="module",
        module="scripts.particles.pic_em_vacuum_wave",
        output_globs=("tst/build/src/pic_em_vacuum_np*-errs.dat",),
        note="EM vacuum proxy convergence from L1 error tables.",
    ),
    PublicationCase(
        case_id="langmuir_frequency_proxy",
        group="entity_core_proxy",
        runner="module",
        module="scripts.particles.pic_langmuir_frequency_proxy",
        output_globs=("tst/build/src/bin/pic_langmuir_freq_proxy_np*.*.bin",),
        note="Frequency extraction from integrated particle-current traces.",
    ),
    PublicationCase(
        case_id="two_stream_growth_proxy",
        group="entity_core_proxy",
        runner="module",
        module="scripts.particles.pic_two_stream_growth_proxy",
        output_globs=("tst/build/src/bin/pic_two_stream_growth_proxy_np*.*.bin",),
        note="Mode-growth envelope from deposited charge density.",
    ),
    PublicationCase(
        case_id="weibel_growth_proxy",
        group="entity_core_proxy",
        runner="module",
        module="scripts.particles.pic_weibel_growth_proxy",
        output_globs=("tst/build/src/bin/pic_weibel_growth_proxy_np*.*.bin",),
        note="Mode-growth envelope from transverse current.",
    ),
    PublicationCase(
        case_id="bell_growth_proxy",
        group="benchmarks_proxy",
        runner="module",
        module="scripts.particles.pic_bell_growth_proxy",
        output_globs=("tst/build/src/bin/pic_bell_proxy_*.*.bin",),
        note="Coupled-vs-uncoupled Bell-like magnetic growth proxy.",
    ),
    PublicationCase(
        case_id="multispecies_backreaction_oscillation",
        group="benchmarks_proxy",
        runner="module",
        module="scripts.particles.pic_multispecies_backreaction_oscillation",
        output_globs=("tst/build/src/bin/pic_mso_*.*.bin",),
        note="Uniform/SMR/AMR-proxy oscillation and energy-drift trends.",
    ),
    PublicationCase(
        case_id="crsi_deltaf_proxy",
        group="benchmarks_proxy",
        runner="module",
        module="scripts.particles.pic_crsi_deltaf_proxy",
        output_globs=("tst/build/src/bin/pic_crsi_deltaf_*.*.bin",),
        note="CRSI polarization growth with delta-f noise suppression.",
    ),
    PublicationCase(
        case_id="crpai_polarization_proxy",
        group="benchmarks_proxy",
        runner="module",
        module="scripts.particles.pic_crpai_polarization_proxy",
        output_globs=("tst/build/src/bin/pic_crpai_*.*.bin",),
        note="CRPAI branch-selection and polarization behavior.",
    ),
    PublicationCase(
        case_id="expanding_box_anisotropy_proxy",
        group="benchmarks_proxy",
        runner="module",
        module="scripts.particles.pic_expanding_box_anisotropy_proxy",
        output_globs=("tst/build/src/bin/pic_box_*.*.bin",),
        note="Expanding/compressing anisotropy trend split.",
    ),
    # ---------------------------------------------------------------------
    # Publication-physics variants (preserve proxy suite, add stronger runs)
    # ---------------------------------------------------------------------
    PublicationCase(
        case_id="entity_deposit_mink_publication",
        group="entity_core_publication",
        runner="module",
        module="scripts.particles.pic_entity_deposit_mink",
        output_globs=("tst/build/src/bin/pic_entity_dep_mink_*.*.bin",),
        note="Publication panel source: Entity parity Minkowski support map.",
    ),
    PublicationCase(
        case_id="entity_deposit_reflect_publication",
        group="entity_core_publication",
        runner="module",
        module="scripts.particles.pic_entity_deposit_reflect",
        output_globs=("tst/build/src/bin/pic_entity_dep_reflect_*.*.bin",),
        note="Publication panel source: Entity parity reflecting support map.",
    ),
    PublicationCase(
        case_id="em_vacuum_wave_publication",
        group="entity_core_publication",
        runner="module",
        module="scripts.particles.pic_em_vacuum_wave_publication",
        output_globs=("tst/build/src/pic_em_vacuum_pub_np*-errs.dat",),
        note="Publication convergence sweep with 5+ resolutions and MPI parity.",
    ),
    PublicationCase(
        case_id="langmuir_frequency_publication",
        group="entity_core_publication",
        runner="module",
        module="scripts.particles.pic_langmuir_frequency_proxy",
        output_globs=("tst/build/src/bin/pic_langmuir_freq_proxy_np*.*.bin",),
        note="Publication panel source: Langmuir oscillation with analytic overlay.",
    ),
    PublicationCase(
        case_id="two_stream_growth_publication",
        group="entity_core_publication",
        runner="module",
        module="scripts.particles.pic_two_stream_growth_publication",
        output_globs=("tst/build/src/bin/pic_two_stream_pub_np*.*.bin",),
        note="Coupled MHD-PIC two-stream growth publication regime.",
    ),
    PublicationCase(
        case_id="weibel_growth_publication",
        group="entity_core_publication",
        runner="module",
        module="scripts.particles.pic_weibel_growth_publication",
        output_globs=("tst/build/src/bin/pic_weibel_pub_np*.*.bin",),
        note="Coupled MHD-PIC Weibel growth publication regime.",
    ),
    PublicationCase(
        case_id="bell_growth_publication",
        group="benchmarks_publication",
        runner="module",
        module="scripts.particles.pic_bell_growth_publication",
        output_globs=("tst/build/src/bin/pic_bell_pub_*.*.bin",),
        note="Extended Bell growth with auto-fit linear window and MPI overlays.",
    ),
    PublicationCase(
        case_id="multispecies_backreaction_publication",
        group="benchmarks_publication",
        runner="module",
        module="scripts.particles.pic_multispecies_backreaction_publication",
        output_globs=("tst/build/src/bin/pic_mso_pub_*.*.bin",),
        note=(
            "Long uniform/SMR multispecies oscillations plus conservative AMR "
            "drift check."
        ),
    ),
    PublicationCase(
        case_id="crsi_deltaf_publication",
        group="benchmarks_publication",
        runner="module",
        module="scripts.particles.pic_crsi_deltaf_publication",
        output_globs=("tst/build/src/bin/pic_crsi_pub_*.*.bin",),
        note="CRSI publication regime with branch growth/rank spread diagnostics.",
    ),
    PublicationCase(
        case_id="crpai_polarization_publication",
        group="benchmarks_publication",
        runner="module",
        module="scripts.particles.pic_crpai_polarization_publication",
        output_globs=("tst/build/src/bin/pic_crpai_pub_*.*.bin",),
        note="CRPAI publication regime with prolate/oblate branch asymmetry checks.",
    ),
    PublicationCase(
        case_id="expanding_box_anisotropy_publication",
        group="benchmarks_publication",
        runner="module",
        module="scripts.particles.pic_expanding_box_anisotropy_proxy",
        output_globs=("tst/build/src/bin/pic_box_*.*.bin",),
        note="Publication panel source: expanding/compressing anisotropy trends.",
    ),
    # ---------------------------------------------------------------------
    # CR-shock MHD-PIC storyline
    # ---------------------------------------------------------------------
    PublicationCase(
        case_id="amr_shock_lb_smoke",
        group="shock_story",
        runner="module",
        module="scripts.particles.pic_amr_shock_lb_smoke",
        output_globs=("tst/build/src/bin/pic_amr_shock_lb_*.*.bin",),
        note="Fast regression smoke anchor for shock/AMR/LB behavior.",
    ),
    PublicationCase(
        case_id="amr_shock_publication_local",
        group="shock_story",
        runner="deck",
        input_deck="inputs/tests/pic_amr_shock_lb_publication_local.athinput",
        output_globs=(
            "tst/build/src/bin/pic_amr_shock_pub_local.*.bin",
            "tst/build/src/pvtk/pic_amr_shock_pub_local.*.part.vtk",
        ),
        ranks_local=2,
        ranks_hpc=8,
        profile="both",
        note="Publication local profile with richer fields plus particle dumps.",
    ),
    PublicationCase(
        case_id="amr_shock_publication_hpc",
        group="shock_story",
        runner="deck",
        input_deck="inputs/tests/pic_amr_shock_lb_publication_hpc.athinput",
        output_globs=(
            "tst/build/src/bin/pic_amr_shock_pub_hpc.*.bin",
            "tst/build/src/pvtk/pic_amr_shock_pub_hpc.*.part.vtk",
        ),
        ranks_local=4,
        ranks_hpc=64,
        profile="hpc",
        note="Heavy production profile intended for GPU-capable clusters.",
    ),
)


def select_cases(groups: Iterable[str], tier: str) -> list[PublicationCase]:
    """Return manifest cases filtered by group and execution tier."""
    gset = {g.strip() for g in groups if g.strip()}
    out = []
    for case in PUBLICATION_CASES:
        if gset and case.group not in gset:
            continue
        if tier == "local" and case.profile == "hpc":
            continue
        if tier == "hpc" and case.profile == "local":
            continue
        out.append(case)
    return out
