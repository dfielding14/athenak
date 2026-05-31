#!/usr/bin/env python3
"""Materialize exact source-local Q-023 Bell preparation deck variants."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import tempfile
from typing import Any, Iterable

from tst.publication import analyze_q023_paper_bell_linear as bell


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = bell.CAMPAIGN_ID
SCHEMA_VERSION = 1
ARTIFACT_ROLE = "source_local_q023_paper_bell_linear_deck_preparation_only"
QUALIFICATION_EFFECT = "none_source_local_preparation_only"
SECTION52_QUALIFICATION = "not_claimed"
FRONTIER_AUTHORIZATION = "not_bound"
AUTHORIZED_OUTPUT_PARENT = REPO_ROOT / "tst/.codex"
DEFAULT_OUTPUT_ROOT = (
    AUTHORIZED_OUTPUT_PARENT / "q023-paper-bell-linear-materialized-variants"
)

EPSILON_VALUES = bell.EPSILON_VALUES
DIMENSIONS = bell.DIMENSIONS
PPC_VALUES = (2, 8, 32)
RESOLUTION_SCALES = (0.5, 1.0, 2.0)
TIMESTEP_SCALES = (0.5, 1.0, 2.0)
EXPECTED_VARIANT_COUNT = (
    len(DIMENSIONS)
    * len(EPSILON_VALUES)
    * len(RESOLUTION_SCALES)
    * len(TIMESTEP_SCALES)
    * len(PPC_VALUES)
)

SOURCE_DECKS = dict(bell.DECKS)
SOURCE_DECK_SHA256 = {
    1: "b84fb1b046d3174c3ce6bdd5f9067173dccd0636557d57c047331bf0a42eb4d5",
    2: "b8a11ce7f7227210392de133bf52f71c265d86ab5ade73bf89a45a6449a3a23c",
    3: "a0fda35a71542e3c411c574a4c5069cbbae3776f02aa10dd44176a3102664103",
}

LOADING_POLICY = {
    "cr_distribution": "center",
    "randomized_seed_materialization": "not_permitted_without_review",
    "review_status": (
        "unresolved_reviewer_boundary_retain_deterministic_centered_loading_"
        "with_deterministic_tolerance_handling_or_approve_randomized_seeded_loading"
    ),
}
TIMESTEP_POLICY = (
    "source_local_cfl_scaling_around_reviewed_base_clean_candidate_timestep_"
    "freeze_still_open"
)
RESOLUTION_POLICY = (
    "scale_active_global_and_meshblock_cell_counts_preserve_bounds_and_"
    "physical_1d_transverse_invariant_nx2_4_carrier"
)

_REQUEST_KEYS = {
    "schema_version",
    "campaign_id",
    "artifact_role",
    "qualification_effect",
    "section52_qualification",
    "qualifying_evidence",
    "frontier_authorization",
    "loading_policy",
    "timestep_policy",
    "resolution_policy",
    "variants",
}
_SELECTOR_KEYS = {
    "variant_id",
    "dimension",
    "epsilon",
    "resolution_scale",
    "timestep_scale",
    "ppc",
}
_BLOCK_HEADER = re.compile(r"^\s*<([^>]+)>\s*(?:#.*)?$")
_PARAMETER_LINE = re.compile(
    r"^(?P<prefix>\s*)(?P<name>[^#=\s]+)(?P<between>\s*=\s*)"
    r"(?P<value>.*?)(?P<suffix>\s*(?:#.*)?)$"
)


class ContractError(ValueError):
    """Raised when Q-023 preparation materialization fails closed."""


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _sha256(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _json_text(document: dict[str, Any]) -> str:
    return json.dumps(
        document, indent=2, sort_keys=True, allow_nan=False
    ) + "\n"


def _finite_float(label: str, value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{label} must be a finite number")
    measured = float(value)
    if not math.isfinite(measured):
        raise ContractError(f"{label} must be a finite number")
    return measured


def _require_close(label: str, measured: float, expected: float) -> None:
    if not math.isclose(measured, expected, rel_tol=1.0e-14, abs_tol=1.0e-14):
        raise ContractError(
            f"{label}: expected {expected!r}, measured {measured!r}"
        )


def _float_token(value: float) -> str:
    return repr(float(value)).replace(".", "p")


def _variant_id(
    dimension: int,
    epsilon: float,
    resolution_scale: float,
    timestep_scale: float,
    ppc: int,
) -> str:
    return (
        f"Q023-SOURCE-LOCAL-PREPARATION-{dimension}D-"
        f"EPSILON-{_float_token(epsilon)}-"
        f"RESOLUTION-{_float_token(resolution_scale)}-"
        f"TIMESTEP-{_float_token(timestep_scale)}-PPC-{ppc}"
    )


def _basename(
    dimension: int,
    epsilon: float,
    resolution_scale: float,
    timestep_scale: float,
    ppc: int,
) -> str:
    return (
        f"pic_q023_paper_bell_linear_{dimension}d_"
        f"epsilon_{_float_token(epsilon)}_"
        f"resolution_{_float_token(resolution_scale)}_"
        f"timestep_{_float_token(timestep_scale)}_ppc_{ppc}_preparation"
    )


def _expected_selectors() -> Iterable[dict[str, Any]]:
    for dimension in DIMENSIONS:
        for epsilon in EPSILON_VALUES:
            for resolution_scale in RESOLUTION_SCALES:
                for timestep_scale in TIMESTEP_SCALES:
                    for ppc in PPC_VALUES:
                        yield {
                            "variant_id": _variant_id(
                                dimension,
                                epsilon,
                                resolution_scale,
                                timestep_scale,
                                ppc,
                            ),
                            "dimension": dimension,
                            "epsilon": epsilon,
                            "resolution_scale": resolution_scale,
                            "timestep_scale": timestep_scale,
                            "ppc": ppc,
                        }


def build_materialization_request() -> dict[str, Any]:
    """Build the complete fixed source-local preparation request."""
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "section52_qualification": SECTION52_QUALIFICATION,
        "qualifying_evidence": False,
        "frontier_authorization": FRONTIER_AUTHORIZATION,
        "loading_policy": dict(LOADING_POLICY),
        "timestep_policy": TIMESTEP_POLICY,
        "resolution_policy": RESOLUTION_POLICY,
        "variants": list(_expected_selectors()),
    }


def _refuse_qualification_claims(document: dict[str, Any]) -> None:
    boundaries = {
        "qualification_effect": QUALIFICATION_EFFECT,
        "section52_qualification": SECTION52_QUALIFICATION,
        "qualifying_evidence": False,
        "frontier_authorization": FRONTIER_AUTHORIZATION,
    }
    for key, expected in boundaries.items():
        if document.get(key) != expected:
            raise ContractError("Q-023 materializer refuses qualification claims")
    for key in document:
        if key not in _REQUEST_KEYS and (
            "qualif" in key.lower() or "authoriz" in key.lower()
        ):
            raise ContractError("Q-023 materializer refuses qualification claims")


def _validated_selector(selector: Any) -> dict[str, Any]:
    if not isinstance(selector, dict):
        raise ContractError("Q-023 materialization variant must be an object")
    for key in selector:
        if key not in _SELECTOR_KEYS and (
            "qualif" in key.lower() or "authoriz" in key.lower()
        ):
            raise ContractError("Q-023 materializer refuses qualification claims")
    if set(selector) != _SELECTOR_KEYS:
        raise ContractError("Q-023 materialization variant fields do not match")

    dimension = selector["dimension"]
    if type(dimension) is not int or dimension not in DIMENSIONS:
        raise ContractError("Q-023 materialization dimension is outside the grid")
    epsilon = _finite_float("epsilon", selector["epsilon"])
    resolution_scale = _finite_float(
        "resolution_scale", selector["resolution_scale"]
    )
    timestep_scale = _finite_float("timestep_scale", selector["timestep_scale"])
    ppc = selector["ppc"]
    if type(ppc) is not int or ppc not in PPC_VALUES:
        raise ContractError("Q-023 materialization PPC is outside the grid")
    if epsilon not in EPSILON_VALUES:
        raise ContractError("Q-023 materialization epsilon is outside the grid")
    if resolution_scale not in RESOLUTION_SCALES:
        raise ContractError("Q-023 materialization resolution is outside the grid")
    if timestep_scale not in TIMESTEP_SCALES:
        raise ContractError("Q-023 materialization timestep is outside the grid")
    expected_id = _variant_id(
        dimension, epsilon, resolution_scale, timestep_scale, ppc
    )
    if selector["variant_id"] != expected_id:
        raise ContractError("Q-023 materialization variant ID mismatch")
    return {
        "variant_id": expected_id,
        "dimension": dimension,
        "epsilon": epsilon,
        "resolution_scale": resolution_scale,
        "timestep_scale": timestep_scale,
        "ppc": ppc,
    }


def validate_materialization_request(request: Any) -> list[dict[str, Any]]:
    """Validate the fixed preparation grid and reject any broadened claim."""
    if not isinstance(request, dict):
        raise ContractError("Q-023 materialization request must be an object")
    _refuse_qualification_claims(request)
    if set(request) != _REQUEST_KEYS:
        raise ContractError("Q-023 materialization request fields do not match")
    if (
        request["schema_version"] != SCHEMA_VERSION
        or request["campaign_id"] != CAMPAIGN_ID
        or request["artifact_role"] != ARTIFACT_ROLE
    ):
        raise ContractError("Q-023 materialization request identity mismatch")
    if request["loading_policy"] != LOADING_POLICY:
        raise ContractError("Q-023 centered-loading review boundary mismatch")
    if request["timestep_policy"] != TIMESTEP_POLICY:
        raise ContractError("Q-023 source-local timestep policy mismatch")
    if request["resolution_policy"] != RESOLUTION_POLICY:
        raise ContractError("Q-023 source-local resolution policy mismatch")
    if not isinstance(request["variants"], list):
        raise ContractError("Q-023 materialization variants must be a list")

    selectors = []
    seen_keys = set()
    seen_ids = set()
    for raw_selector in request["variants"]:
        selector = _validated_selector(raw_selector)
        key = (
            selector["dimension"],
            selector["epsilon"],
            selector["resolution_scale"],
            selector["timestep_scale"],
            selector["ppc"],
        )
        if key in seen_keys or selector["variant_id"] in seen_ids:
            raise ContractError("duplicate Q-023 materialization variant")
        seen_keys.add(key)
        seen_ids.add(selector["variant_id"])
        selectors.append(selector)

    expected = list(_expected_selectors())
    if selectors != expected:
        raise ContractError(
            "Q-023 materialization grid is incomplete, reordered, or contains extras"
        )
    return selectors


def _validated_source_decks() -> dict[int, dict[str, Any]]:
    analyzed = {
        deck["dimension"]: deck for deck in bell.validate_source_local_candidate_decks()
    }
    decks = {}
    for dimension in DIMENSIONS:
        path = SOURCE_DECKS[dimension]
        digest = _sha256(path)
        if digest != SOURCE_DECK_SHA256[dimension]:
            raise ContractError(
                f"Q-023 reviewed {dimension}D source deck digest mismatch"
            )
        blocks = bell.parse_athinput(path)
        if blocks["particles"]["cr_distribution"] != "center":
            raise ContractError("Q-023 reviewed source deck is not centered loading")
        decks[dimension] = {
            "path": path,
            "relative_path": str(path.relative_to(REPO_ROOT)),
            "sha256": digest,
            "text": path.read_text(encoding="utf-8"),
            "blocks": blocks,
            "analyzed": analyzed[dimension],
        }
    return decks


def _render_float(value: float, original: str | None = None) -> str:
    measured = float(value)
    if original is not None:
        try:
            if math.isclose(
                measured, float(original), rel_tol=1.0e-14, abs_tol=1.0e-14
            ):
                return original
        except ValueError:
            pass
    return repr(measured)


def _scaled_count(original: str, scale: float) -> str:
    value = int(original)
    scaled = value * scale
    if not scaled.is_integer() or scaled < 1:
        raise ContractError("Q-023 scaled resolution does not produce cell counts")
    return str(int(scaled))


def _rewrite_deck(
    source_text: str, replacements: dict[tuple[str, str], str]
) -> str:
    rendered = []
    replaced = set()
    current_block = None
    for raw_line in source_text.splitlines(keepends=True):
        content = raw_line[:-1] if raw_line.endswith("\n") else raw_line
        newline = "\n" if raw_line.endswith("\n") else ""
        header = _BLOCK_HEADER.fullmatch(content)
        if header is not None:
            current_block = header.group(1).strip()
            rendered.append(raw_line)
            continue
        parameter = _PARAMETER_LINE.fullmatch(content)
        if parameter is None or current_block is None:
            rendered.append(raw_line)
            continue
        key = (current_block, parameter.group("name"))
        if key not in replacements:
            rendered.append(raw_line)
            continue
        if key in replaced:
            raise ContractError(f"Q-023 source deck contains duplicate {key!r}")
        replaced.add(key)
        rendered.append(
            parameter.group("prefix")
            + parameter.group("name")
            + parameter.group("between")
            + replacements[key]
            + parameter.group("suffix")
            + newline
        )
    if replaced != set(replacements):
        missing = sorted(set(replacements) - replaced)
        raise ContractError(
            f"Q-023 source deck replacement targets are missing: {missing}"
        )
    return "".join(rendered)


def _active_axes(dimension: int) -> tuple[int, ...]:
    return (1,) if dimension == 1 else tuple(range(1, dimension + 1))


def _render_variant(
    selector: dict[str, Any], source: dict[str, Any]
) -> tuple[dict[str, Any], str]:
    blocks = source["blocks"]
    particles = blocks["particles"]
    time = blocks["time"]
    basis = bell._mode_basis(selector["dimension"])[0]
    stream_speed = 1.0 / selector["epsilon"]
    stream_velocity = [stream_speed * float(component) for component in basis]
    light_speed = float(
        blocks["q023_paper_bell_linear"]["c_over_v_cr"]
    ) * stream_speed
    base_ppc = float(particles["ppc"])
    base_qscale = float(particles["deposit_qscale"])
    deposit_qscale = base_ppc * base_qscale / selector["ppc"]
    cfl_number = float(time["cfl_number"]) * selector["timestep_scale"]
    basename = _basename(
        selector["dimension"],
        selector["epsilon"],
        selector["resolution_scale"],
        selector["timestep_scale"],
        selector["ppc"],
    )

    replacements = {
        ("job", "basename"): basename,
        ("time", "cfl_number"): _render_float(
            cfl_number, time["cfl_number"]
        ),
        ("particles", "ppc"): _render_float(
            selector["ppc"], particles["ppc"]
        ),
        ("particles", "deposit_qscale"): _render_float(
            deposit_qscale, particles["deposit_qscale"]
        ),
        ("particles", "cr_vx0"): _render_float(
            stream_velocity[0], particles["cr_vx0"]
        ),
        ("particles", "cr_vy0"): _render_float(
            stream_velocity[1], particles["cr_vy0"]
        ),
        ("particles", "cr_vz0"): _render_float(
            stream_velocity[2], particles["cr_vz0"]
        ),
        ("particles", "pic_cr_light_speed"): _render_float(
            light_speed, particles["pic_cr_light_speed"]
        ),
        ("q023_paper_bell_linear", "epsilon"): _render_float(
            selector["epsilon"],
            blocks["q023_paper_bell_linear"]["epsilon"],
        ),
    }
    for block_name in ("mesh", "meshblock"):
        for axis in _active_axes(selector["dimension"]):
            field = f"nx{axis}"
            replacements[(block_name, field)] = _scaled_count(
                blocks[block_name][field], selector["resolution_scale"]
            )

    deck_text = _rewrite_deck(source["text"], replacements)
    deck_bytes = deck_text.encode("utf-8")
    deck_path = f"decks/{basename}.athinput"
    record = {
        **selector,
        "source_deck_path": source["relative_path"],
        "source_deck_sha256": source["sha256"],
        "deck_path": deck_path,
        "deck_sha256": _sha256_bytes(deck_bytes),
        "basename": basename,
        "stream_velocity": stream_velocity,
        "pic_cr_light_speed": light_speed,
        "deposit_qscale": deposit_qscale,
        "cfl_number": cfl_number,
        "cr_distribution": "center",
        "qualification_effect": QUALIFICATION_EFFECT,
        "section52_qualification": SECTION52_QUALIFICATION,
        "qualifying_evidence": False,
    }
    _validate_rendered_deck(record, deck_text)
    return record, deck_text


def _validate_rendered_deck(record: dict[str, Any], deck_text: str) -> None:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "variant.athinput"
        path.write_text(deck_text, encoding="utf-8")
        blocks = bell.parse_athinput(path)
    particles = blocks["particles"]
    metadata = blocks["q023_paper_bell_linear"]
    if particles["cr_distribution"] != "center":
        raise ContractError("Q-023 rendered deck changed centered loading")
    if metadata["timestep"] != "open_clean_candidate_timestep_freeze":
        raise ContractError("Q-023 rendered deck claimed a frozen timestep")
    if metadata["epsilon"] != _render_float(
        record["epsilon"], metadata["epsilon"]
    ):
        raise ContractError("Q-023 rendered epsilon is inconsistent")

    stream = [
        float(particles["cr_vx0"]),
        float(particles["cr_vy0"]),
        float(particles["cr_vz0"]),
    ]
    speed = math.sqrt(sum(value * value for value in stream))
    ppc = float(particles["ppc"])
    qscale = float(particles["deposit_qscale"])
    charge = float(blocks["species0"]["charge"])
    light_speed = float(particles["pic_cr_light_speed"])
    b_g = float(metadata["b_g"])
    k0 = float(metadata["k0"])
    _require_close("Q-023 rendered epsilon physics", speed, 1.0 / record["epsilon"])
    _require_close("Q-023 rendered PPC", ppc, record["ppc"])
    _require_close("Q-023 rendered deposit_qscale", qscale, record["deposit_qscale"])
    _require_close(
        "Q-023 rendered artificial light speed",
        light_speed,
        1000.0 * speed,
    )
    _require_close(
        "Q-023 rendered j_CR",
        ppc * qscale * charge * speed,
        2.0 * b_g * light_speed * k0,
    )


def render_variant_deck(selector: Any) -> tuple[dict[str, Any], str]:
    """Render one admitted source-local preparation selector."""
    validated = _validated_selector(selector)
    source = _validated_source_decks()[validated["dimension"]]
    return _render_variant(validated, source)


def build_materialization_manifest(
    request: Any,
) -> tuple[dict[str, Any], dict[str, str], str]:
    """Build deterministic manifest bytes and exact rendered deck text."""
    selectors = validate_materialization_request(request)
    request_text = _json_text(request)
    sources = _validated_source_decks()
    variants = []
    decks = {}
    for selector in selectors:
        record, deck_text = _render_variant(
            selector, sources[selector["dimension"]]
        )
        if record["deck_path"] in decks:
            raise ContractError("duplicate Q-023 materialized deck path")
        decks[record["deck_path"]] = deck_text
        variants.append(record)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "section52_qualification": SECTION52_QUALIFICATION,
        "qualifying_evidence": False,
        "frontier_authorization": FRONTIER_AUTHORIZATION,
        "loading_policy": dict(LOADING_POLICY),
        "timestep_policy": TIMESTEP_POLICY,
        "resolution_policy": RESOLUTION_POLICY,
        "materialization_request_path": "materialization_request.json",
        "materialization_request_sha256": _sha256_bytes(
            request_text.encode("utf-8")
        ),
        "source_decks": [
            {
                "dimension": dimension,
                "path": sources[dimension]["relative_path"],
                "sha256": sources[dimension]["sha256"],
            }
            for dimension in DIMENSIONS
        ],
        "variant_count": len(variants),
        "variants": variants,
        "status": "source_local_preparation_only_no_qualification_claim",
    }
    return manifest, decks, request_text


def _safe_output_root(output_root: Path) -> Path:
    root = Path(output_root)
    if not root.is_absolute():
        raise ContractError("Q-023 materialization output root must be absolute")
    try:
        authorized_parent = AUTHORIZED_OUTPUT_PARENT.resolve(strict=True)
        parent = root.parent.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ContractError(
            "Q-023 materialization output root parent must already exist"
        ) from error
    if not authorized_parent.is_dir() or not parent.is_dir():
        raise ContractError(
            "Q-023 materialization output root parent must be a directory"
        )
    try:
        root.relative_to(authorized_parent)
    except ValueError as error:
        raise ContractError(
            "Q-023 materialization output root is outside tst/.codex"
        ) from error
    if parent != authorized_parent or root.parent != authorized_parent:
        raise ContractError(
            "Q-023 materialization output root must be a direct tst/.codex child"
        )
    if os.path.lexists(root):
        raise ContractError("Q-023 materialization output root already exists")
    return root


def _write_new_text(parent_fd: int, name: str, text: str) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    flags |= getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(name, flags, 0o600, dir_fd=parent_fd)
    try:
        content = text.encode("utf-8")
        written = 0
        while written < len(content):
            written += os.write(fd, content[written:])
    finally:
        os.close(fd)


def _remove_tree_at(parent_fd: int, name: str) -> None:
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    child_fd = os.open(name, flags, dir_fd=parent_fd)
    try:
        for child in os.listdir(child_fd):
            mode = os.stat(
                child, dir_fd=child_fd, follow_symlinks=False
            ).st_mode
            if stat.S_ISDIR(mode):
                _remove_tree_at(child_fd, child)
            else:
                os.unlink(child, dir_fd=child_fd)
    finally:
        os.close(child_fd)
    os.rmdir(name, dir_fd=parent_fd)


def materialize_variant_decks(
    output_root: Path, request: Any | None = None
) -> dict[str, Any]:
    """Write exact preparation decks and manifests below a safe new root."""
    root = _safe_output_root(output_root)
    if request is None:
        request = build_materialization_request()
    manifest, decks, request_text = build_materialization_manifest(request)
    manifest_text = _json_text(manifest)
    parent_flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    parent_fd = os.open(root.parent, parent_flags)
    reserved = False
    try:
        try:
            os.mkdir(root.name, 0o700, dir_fd=parent_fd)
        except FileExistsError as error:
            raise ContractError(
                "Q-023 materialization output root already exists"
            ) from error
        reserved = True
        output_fd = os.open(root.name, parent_flags, dir_fd=parent_fd)
        try:
            os.mkdir("decks", 0o700, dir_fd=output_fd)
            decks_fd = os.open("decks", parent_flags, dir_fd=output_fd)
            try:
                _write_new_text(
                    output_fd, "materialization_request.json", request_text
                )
                for relative, deck_text in decks.items():
                    prefix, name = relative.split("/", 1)
                    if prefix != "decks" or "/" in name:
                        raise ContractError(
                            "Q-023 materialized deck path is not normalized"
                        )
                    _write_new_text(decks_fd, name, deck_text)
            finally:
                os.close(decks_fd)
            _write_new_text(
                output_fd, "materialization_manifest.json", manifest_text
            )
        finally:
            os.close(output_fd)
        reserved = False
    except BaseException:
        if reserved:
            _remove_tree_at(parent_fd, root.name)
        raise
    finally:
        os.close(parent_fd)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", type=Path)
    parser.add_argument("--print-request", action="store_true")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    if args.print_request:
        if args.request is not None:
            parser.error("--print-request does not accept --request")
        print(_json_text(build_materialization_request()), end="")
        return
    request = (
        build_materialization_request()
        if args.request is None
        else json.loads(args.request.read_text(encoding="utf-8"))
    )
    result = materialize_variant_decks(args.output_root, request)
    print(_json_text(result), end="")


if __name__ == "__main__":
    main()
