import json
from pathlib import Path
import sys

CONTROL_PLANE_DIR = Path(__file__).resolve().parent / "frontier_control_plane"
sys.path.insert(0, str(CONTROL_PLANE_DIR))

import generate_prepared_pic_artifact_inventory as generator


def test_q019_campaign_is_in_prepared_artifact_inventory() -> None:
    inventory = generator.prepared_artifact_inventory()
    deck_paths = {record["path"] for record in inventory["paper_decks"]}
    analyzer_paths = {record["path"] for record in inventory["analyzers"]}

    manifest_path = generator.Q019_DECK_ROOT / "deck_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_decks = {record["path"] for record in manifest["decks"]}
    controller_manifest = json.loads(
        (
            generator.Q019_RUNTIME_CONTROLLER_DECK_ROOT / "deck_manifest.json"
        ).read_text(encoding="utf-8")
    )
    expected_controller_decks = {
        (
            generator.Q019_RUNTIME_CONTROLLER_DECK_ROOT
            / str(record["filename"])
        ).relative_to(generator.REPO_ROOT).as_posix()
        for record in controller_manifest["artifacts"]
    }

    assert len(expected_decks) == 77
    assert expected_decks <= deck_paths
    assert len(expected_controller_decks) == 25
    assert expected_controller_decks <= deck_paths
    assert (
        "tst/publication/analyze_q019_physics_first_nonlinear_bell_successor_v2.py"
        in analyzer_paths
    )


def test_committed_inventory_matches_generator() -> None:
    committed = json.loads(generator.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
    assert committed == generator.prepared_artifact_inventory()
