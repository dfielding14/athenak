"""Manifest for exploratory PIC/MHD-PIC engineering artifacts.

Historical identifiers containing ``publication`` are retained for CLI and
archive compatibility only. Every checked-in case is explicitly unqualified:
none is Sun and Bai (2023) reproduction evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
from typing import Tuple


ENGINEERING_PROXY = "engineering_proxy"

LEGACY_GROUP_ALIASES = {
    "entity_core_proxy": "entity_core_engineering",
    "benchmarks_proxy": "benchmark_engineering",
    "entity_core_publication": "extended_entity_engineering",
    "benchmarks_publication": "extended_benchmark_engineering",
    "shock_story": "shock_engineering",
}


@dataclass(frozen=True)
class ArtifactCase:
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
    evidence_class: str = ENGINEERING_PROXY
    note: str = ""

    @property
    def legacy_identifier(self) -> bool:
        """Return whether the compatibility identifier predates reclassification."""
        return "publication" in self.case_id

    @property
    def physical_model(self) -> str:
        """Return a scoped identifier for the unqualified model exercised."""
        cid = self.case_id
        if "entity_deposit" in cid:
            return "entity_inspired_deposition_parity"
        if "em_vacuum" in cid:
            return "mhd_linear_wave_with_inactive_particles"
        if "langmuir" in cid:
            return "nonrelativistic_uniform_b_gyrofrequency_anchor"
        if "two_stream" in cid:
            return "two_stream_like_engineering_control"
        if "weibel" in cid:
            return "transverse_current_weibel_like_engineering_control"
        if "bell" in cid:
            return "bell_like_current_coupling_engineering_control"
        if "multispecies" in cid:
            return "multispecies_backreaction_oscillation_engineering_control"
        if "crsi" in cid:
            return "crsi_quiet_start_engineering_control"
        if "crpai" in cid:
            return "crpai_branch_selectivity_engineering_control"
        if "expanding_box" in cid:
            return "particle_only_box_scaling_engineering_control"
        if "shock" in cid:
            return "orszag_tang_shock_rich_engineering_stress"
        return "exploratory_unqualified"


ARTIFACT_CASES: Tuple[ArtifactCase, ...] = (
    # ---------------------------------------------------------------------
    # Proxy-regression anchors (unchanged physics and thresholds)
    # ---------------------------------------------------------------------
    ArtifactCase(
        case_id="entity_deposit_mink",
        group="entity_core_engineering",
        runner="module",
        module="scripts.particles.pic_entity_deposit_mink",
        output_globs=("tst/build/src/bin/pic_entity_dep_mink_*.*.bin",),
        note="Entity parity: Minkowski deposit invariance.",
    ),
    ArtifactCase(
        case_id="entity_deposit_reflect",
        group="entity_core_engineering",
        runner="module",
        module="scripts.particles.pic_entity_deposit_reflect",
        output_globs=("tst/build/src/bin/pic_entity_dep_reflect_*.*.bin",),
        note="Entity parity: reflecting-boundary deposit invariance.",
    ),
    ArtifactCase(
        case_id="em_vacuum_wave",
        group="entity_core_engineering",
        runner="module",
        module="scripts.particles.pic_em_vacuum_wave",
        output_globs=("tst/build/src/pic_em_vacuum_np*-errs.dat",),
        note="EM vacuum proxy convergence from L1 error tables.",
    ),
    ArtifactCase(
        case_id="langmuir_frequency_proxy",
        group="entity_core_engineering",
        runner="module",
        module="scripts.particles.pic_langmuir_frequency_proxy",
        output_globs=("tst/build/src/bin/pic_langmuir_freq_proxy_np*.*.bin",),
        note="Frequency extraction from integrated particle-current traces.",
    ),
    ArtifactCase(
        case_id="two_stream_growth_proxy",
        group="entity_core_engineering",
        runner="module",
        module="scripts.particles.pic_two_stream_growth_proxy",
        output_globs=("tst/build/src/bin/pic_two_stream_growth_proxy_np*.*.bin",),
        note="Mode-growth envelope from deposited charge density.",
    ),
    ArtifactCase(
        case_id="weibel_growth_proxy",
        group="entity_core_engineering",
        runner="module",
        module="scripts.particles.pic_weibel_growth_proxy",
        output_globs=("tst/build/src/bin/pic_weibel_growth_proxy_np*.*.bin",),
        note="Mode-growth envelope from transverse current.",
    ),
    ArtifactCase(
        case_id="bell_growth_proxy",
        group="benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_bell_growth_proxy",
        output_globs=("tst/build/src/bin/pic_bell_proxy_*.*.bin",),
        note="Coupled-vs-uncoupled Bell-like magnetic growth proxy.",
    ),
    ArtifactCase(
        case_id="multispecies_backreaction_oscillation",
        group="benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_multispecies_backreaction_oscillation",
        output_globs=("tst/build/src/bin/pic_mso_*.*.bin",),
        note="Uniform/SMR/AMR-proxy oscillation and energy-drift trends.",
    ),
    ArtifactCase(
        case_id="crsi_deltaf_proxy",
        group="benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_crsi_deltaf_proxy",
        output_globs=("tst/build/src/bin/pic_crsi_deltaf_*.*.bin",),
        note="CRSI polarization-growth quiet-start engineering proxy.",
    ),
    ArtifactCase(
        case_id="crpai_polarization_proxy",
        group="benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_crpai_polarization_proxy",
        output_globs=("tst/build/src/bin/pic_crpai_*.*.bin",),
        note="CRPAI branch-selection and polarization behavior.",
    ),
    ArtifactCase(
        case_id="expanding_box_anisotropy_proxy",
        group="benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_expanding_box_anisotropy_proxy",
        output_globs=("tst/build/src/bin/pic_box_*.*.bin",),
        note="Expanding/compressing anisotropy trend split.",
    ),
    # ---------------------------------------------------------------------
    # Extended engineering variants (preserve proxy suite, add stronger runs)
    # ---------------------------------------------------------------------
    ArtifactCase(
        case_id="entity_deposit_mink_publication",
        group="extended_entity_engineering",
        runner="module",
        module="scripts.particles.pic_entity_deposit_mink",
        output_globs=("tst/build/src/bin/pic_entity_dep_mink_*.*.bin",),
        note="Legacy extended engineering panel: Entity parity Minkowski support map.",
    ),
    ArtifactCase(
        case_id="entity_deposit_reflect_publication",
        group="extended_entity_engineering",
        runner="module",
        module="scripts.particles.pic_entity_deposit_reflect",
        output_globs=("tst/build/src/bin/pic_entity_dep_reflect_*.*.bin",),
        note="Legacy extended engineering panel: Entity parity reflecting support map.",
    ),
    ArtifactCase(
        case_id="em_vacuum_wave_publication",
        group="extended_entity_engineering",
        runner="module",
        module="scripts.particles.pic_em_vacuum_wave_publication",
        output_globs=("tst/build/src/pic_em_vacuum_pub_np*-errs.dat",),
        note=(
            "Legacy extended engineering convergence sweep with 5+ resolutions "
            "and MPI parity."
        ),
    ),
    ArtifactCase(
        case_id="langmuir_frequency_publication",
        group="extended_entity_engineering",
        runner="module",
        module="scripts.particles.pic_langmuir_frequency_proxy",
        output_globs=("tst/build/src/bin/pic_langmuir_freq_proxy_np*.*.bin",),
        note=(
            "Legacy extended engineering panel: nonrelativistic uniform-B "
            "frequency anchor."
        ),
    ),
    ArtifactCase(
        case_id="two_stream_growth_publication",
        group="extended_entity_engineering",
        runner="module",
        module="scripts.particles.pic_two_stream_growth_publication",
        output_globs=("tst/build/src/bin/pic_two_stream_pub_np*.*.bin",),
        note="Legacy extended engineering two-stream-like growth regime.",
    ),
    ArtifactCase(
        case_id="weibel_growth_publication",
        group="extended_entity_engineering",
        runner="module",
        module="scripts.particles.pic_weibel_growth_publication",
        output_globs=("tst/build/src/bin/pic_weibel_pub_np*.*.bin",),
        note="Legacy extended engineering Weibel-like growth regime.",
    ),
    ArtifactCase(
        case_id="bell_growth_publication",
        group="extended_benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_bell_growth_publication",
        output_globs=("tst/build/src/bin/pic_bell_pub_*.*.bin",),
        note="Legacy extended engineering Bell-like growth with auto-fit window.",
    ),
    ArtifactCase(
        case_id="multispecies_backreaction_publication",
        group="extended_benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_multispecies_backreaction_publication",
        output_globs=("tst/build/src/bin/pic_mso_pub_*.*.bin",),
        note=(
            "Long uniform/SMR multispecies oscillations plus conservative AMR "
            "drift check."
        ),
    ),
    ArtifactCase(
        case_id="crsi_deltaf_publication",
        group="extended_benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_crsi_deltaf_publication",
        output_globs=("tst/build/src/bin/pic_crsi_pub_*.*.bin",),
        note="Legacy extended engineering CRSI-like branch-growth diagnostics.",
    ),
    ArtifactCase(
        case_id="crpai_polarization_publication",
        group="extended_benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_crpai_polarization_publication",
        output_globs=("tst/build/src/bin/pic_crpai_pub_*.*.bin",),
        note="Legacy extended engineering CRPAI-like branch-asymmetry checks.",
    ),
    ArtifactCase(
        case_id="expanding_box_anisotropy_publication",
        group="extended_benchmark_engineering",
        runner="module",
        module="scripts.particles.pic_expanding_box_anisotropy_proxy",
        output_globs=("tst/build/src/bin/pic_box_*.*.bin",),
        note="Legacy extended engineering expanding/compressing anisotropy trends.",
    ),
    # ---------------------------------------------------------------------
    # Shock-rich engineering stress storyline
    # ---------------------------------------------------------------------
    ArtifactCase(
        case_id="amr_shock_lb_smoke",
        group="shock_engineering",
        runner="module",
        module="scripts.particles.pic_amr_shock_lb_smoke",
        output_globs=("tst/build/src/bin/pic_amr_shock_lb_*.*.bin",),
        note="Fast regression smoke anchor for shock/AMR/LB behavior.",
    ),
    ArtifactCase(
        case_id="amr_shock_publication_local",
        group="shock_engineering",
        runner="deck",
        input_deck="inputs/tests/pic_amr_shock_lb_publication_local.athinput",
        output_globs=(
            "tst/build/src/bin/pic_amr_shock_pub_local.*.bin",
            "tst/build/src/pvtk/pic_amr_shock_pub_local.*.part.vtk",
        ),
        ranks_local=2,
        ranks_hpc=8,
        profile="both",
        note="Legacy local engineering-stress profile with richer diagnostics.",
    ),
    ArtifactCase(
        case_id="amr_shock_publication_hpc",
        group="shock_engineering",
        runner="deck",
        input_deck="inputs/tests/pic_amr_shock_lb_publication_hpc.athinput",
        output_globs=(
            "tst/build/src/bin/pic_amr_shock_pub_hpc.*.bin",
            "tst/build/src/pvtk/pic_amr_shock_pub_hpc.*.part.vtk",
        ),
        ranks_local=4,
        ranks_hpc=64,
        profile="hpc",
        note="Legacy heavy engineering-stress profile; not production-qualified.",
    ),
)


def canonicalize_groups(groups: Iterable[str]) -> set[str]:
    """Map deprecated group aliases to canonical engineering group names."""
    return {
        LEGACY_GROUP_ALIASES.get(g.strip(), g.strip())
        for g in groups
        if g.strip()
    }


def select_cases(groups: Iterable[str], tier: str) -> list[ArtifactCase]:
    """Return manifest cases filtered by group and execution tier."""
    gset = canonicalize_groups(groups)
    out = []
    for case in ARTIFACT_CASES:
        if gset and case.group not in gset:
            continue
        if tier == "local" and case.profile == "hpc":
            continue
        if tier == "hpc" and case.profile == "local":
            continue
        out.append(case)
    return out
