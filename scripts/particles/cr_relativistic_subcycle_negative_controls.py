#!/usr/bin/env python3
"""Prove that Phase-5 source-contract guards reject weakened implementations."""

from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import cr_relativistic_subcycle_runtime_inspect as phase5


FILES = (
    "src/particles/particles_pushers.cpp",
    "src/particles/particles_tasks.cpp",
    "src/particles/particles.cpp",
    "src/mesh/mesh_refinement.cpp",
    "src/mesh/mesh.cpp",
    "src/pgen/unit_tests/cr_relativistic_coupled_runtime_test.cpp",
)
MUTATIONS = (
    ("ceil_exact_ratio", FILES[0],
     "if (static_cast<Real>(n) < ratio) {++n;}",
     "if (static_cast<Real>(n) <= ratio) {++n;}"),
    ("drop_checked_multiplication", FILES[0],
     "relativistic::CheckedProduct(envelope.alpha_max, envelope.e_max, alpha_e)",
     "true"),
    ("drop_gamma_floor", FILES[0],
     "envelope.w_min - kick", "envelope.w_min"),
    ("drop_midpoint_regather", FILES[0],
     "              RunRelativisticMHDIdealStep(",
     "              RunRelativisticMHDIdealStep_DISABLED("),
    ("drop_post_integrator_refresh", FILES[1],
     "RefreshRelativisticTimestepBound();",
     "RefreshRelativisticTimestepBound_DISABLED();"),
    ("drop_cap_preflight", FILES[0],
     "if (nsub > subcycle_max_steps) {\n    FatalRelativisticBound(",
     "if (false && nsub > subcycle_max_steps) {\n    FatalRelativisticBound("),
    ("refresh_inside_push", FILES[0],
     'MarkRelativisticTimestepBoundDirty();\n  CheckMotionBounds("particle push");',
     'MarkRelativisticTimestepBoundDirty();\n'
     '  RefreshRelativisticTimestepBound();\n'
     '  CheckMotionBounds("particle push");'),
    ("refresh_during_initialization", FILES[5],
     "pmy_mesh_->pmb_pack->ppart->MarkRelativisticTimestepBoundDirty();",
     "pmy_mesh_->pmb_pack->ppart->MarkRelativisticTimestepBoundDirty();\n"
     "  pmy_mesh_->pmb_pack->ppart->RefreshRelativisticTimestepBound();"),
)


def _copy_contract(root, target):
    for name in FILES:
        destination = target / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / name, destination)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--metrics", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    rejected = {}
    with tempfile.TemporaryDirectory(prefix="cr-rel-phase5-negative-") as tmp:
        scratch = Path(tmp)
        for name, filename, old, new in MUTATIONS:
            shutil.rmtree(scratch, ignore_errors=True)
            scratch.mkdir()
            _copy_contract(root, scratch)
            path = scratch / filename
            text = path.read_text()
            if old not in text:
                raise RuntimeError(f"{name}: mutation target missing")
            path.write_text(text.replace(old, new, 1))
            try:
                phase5._assert_source_order(scratch)
            except Exception as exc:  # Each deliberate weakening must be rejected.
                rejected[name] = str(exc)
            else:
                raise RuntimeError(f"{name}: weakened implementation was accepted")
    metrics = {"rejected_mutation_count": len(rejected), "rejections": rejected}
    if args.metrics:
        args.metrics.parent.mkdir(parents=True, exist_ok=True)
        args.metrics.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print(f"CR relativistic Phase-5 mutation controls passed: {len(rejected)} rejected")


if __name__ == "__main__":
    main()
