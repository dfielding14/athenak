#!/opt/cray/pe/python/3.11.7/bin/python3
"""Install an immutable, checksummed Frontier PIC control-plane version."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import uuid

from control_plane_common import CONTROL_PLANE_FILES, inventory_digest
from control_plane_common import make_tree_read_only, remove_tree, sha256
from control_plane_common import require_no_symlink_components_below


SCRIPT_DIR = Path(__file__).resolve().parent


def install(pic_root: Path) -> Path:
    pic_root = pic_root.resolve()
    records = [
        {"path": name, "sha256": sha256(SCRIPT_DIR / name)}
        for name in CONTROL_PLANE_FILES
    ]
    digest = inventory_digest(records)
    destination = pic_root / "control_plane" / digest
    require_no_symlink_components_below(destination, pic_root)
    if destination.exists():
        raise ValueError(f"Control-plane version already exists: {destination}")
    parent = destination.parent
    parent.mkdir(parents=True, exist_ok=True)
    require_no_symlink_components_below(destination, pic_root)
    temporary = parent / f".tmp-{digest}-{uuid.uuid4()}"
    temporary.mkdir()
    try:
        for record in records:
            source = SCRIPT_DIR / str(record["path"])
            shutil.copy2(source, temporary / source.name)
        inventory = {"schema_version": 1, "version": digest, "files": records}
        (temporary / "inventory.json").write_text(
            json.dumps(inventory, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        make_tree_read_only(
            temporary,
            executable_names={
                path.name for path in temporary.iterdir()
                if path.suffix in {".py", ".sh"}
            },
        )
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            remove_tree(temporary)
    print(destination)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pic-root", required=True, type=Path)
    args = parser.parse_args()
    install(args.pic_root.resolve())


if __name__ == "__main__":
    main()
