import io
import json
from pathlib import Path
import sys
import tarfile

CONTROL_PLANE_DIR = Path(__file__).resolve().parent / "frontier_control_plane"
sys.path.insert(0, str(CONTROL_PLANE_DIR))

from control_plane_common import prepared_artifact_manifest_from_source_archive
import generate_prepared_pic_artifact_inventory as generator


def test_committed_inventory_matches_generator() -> None:
    committed = json.loads(generator.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
    assert committed == generator.prepared_artifact_inventory()


def test_committed_inventory_rederives_from_exact_worktree_bytes() -> None:
    inventory = json.loads(generator.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
    paths = [
        generator.DEFAULT_OUTPUT.relative_to(generator.REPO_ROOT).as_posix(),
        *(record["path"] for record in inventory["paper_decks"]),
        *(record["path"] for record in inventory["analyzers"]),
    ]
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode="w") as stream:
        for relative in sorted(set(paths)):
            payload = (generator.REPO_ROOT / relative).read_bytes()
            member = tarfile.TarInfo(relative)
            member.size = len(payload)
            stream.addfile(member, io.BytesIO(payload))

    manifest = prepared_artifact_manifest_from_source_archive(
        archive.getvalue(),
        inventory_path=(
            "tst/publication/frontier_control_plane/"
            "prepared_pic_artifact_inventory.json"
        ),
    )
    assert manifest["paper_decks"] == inventory["paper_decks"]
    assert manifest["analyzers"] == inventory["analyzers"]
