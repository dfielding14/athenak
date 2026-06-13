#!/usr/bin/env python3
"""Generate the committed clean-candidate PIC deck and analyzer inventory."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

if __package__:
    from .control_plane_common import PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
else:
    from control_plane_common import PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
DEFAULT_OUTPUT = SCRIPT_DIR / "prepared_pic_artifact_inventory.json"
Q043_DECK_ROOT = (
    REPO_ROOT / "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle"
)
Q023_DECK_ROOT = (
    REPO_ROOT / "inputs/tests/q023_paper_bell_linear_joverc_predecessor"
)
Q019_DECK_ROOT = (
    REPO_ROOT / "inputs/publication/q019_physics_first_nonlinear_bell_successor_v2"
)
Q019_RUNTIME_CONTROLLER_DECK_ROOT = (
    REPO_ROOT / "inputs/publication/q019_nonlinear_bell_runtime_controller_v1"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(path: Path) -> dict[str, str]:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": _sha256(path),
    }


def prepared_artifact_inventory() -> dict[str, object]:
    paper_decks = sorted(
        [
            *(REPO_ROOT / "inputs/tests").glob("pic*.athinput"),
            *Q043_DECK_ROOT.glob("*.athinput"),
            *Q023_DECK_ROOT.glob("*.athinput"),
            *Q019_DECK_ROOT.glob("*.athinput"),
            *Q019_RUNTIME_CONTROLLER_DECK_ROOT.glob("*.athinput"),
            *(
                REPO_ROOT / path
                for path in PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
            ),
        ]
    )
    analyzers = sorted((REPO_ROOT / "tst/publication").glob("analyze_*.py"))
    if not paper_decks:
        raise ValueError("No PIC input decks were found")
    if not analyzers:
        raise ValueError("No publication analyzers were found")
    return {
        "schema_version": 1,
        "paper_decks": [_record(path) for path in paper_decks],
        "analyzers": [_record(path) for path in analyzers],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    output.relative_to(REPO_ROOT)
    payload = json.dumps(prepared_artifact_inventory(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload, encoding="utf-8")
    print(output.relative_to(REPO_ROOT).as_posix())


if __name__ == "__main__":
    main()
