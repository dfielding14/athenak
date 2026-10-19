#!/usr/bin/env python3
"""Rewrite staged plot-output paths in published metadata to the final root."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
from pathlib import Path
from typing import Any, Iterable


def _write_json_atomic(path: Path, value: Any) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def repath_value(value: Any, campaign_root: Path, final_root: Path) -> Any:
    if isinstance(value, str):
        path = Path(value)
        try:
            relative = path.relative_to(campaign_root)
        except ValueError:
            return value
        parts = relative.parts
        if not parts or not parts[0].isdigit():
            return value
        suffix = Path(*parts[1:]) if len(parts) > 1 else Path()
        return str(final_root / suffix)
    if isinstance(value, list):
        return [repath_value(item, campaign_root, final_root) for item in value]
    if isinstance(value, dict):
        return {
            key: repath_value(item, campaign_root, final_root)
            for key, item in value.items()
        }
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata_dir", type=Path)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--final-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    return parser


def _rewrite_one(task: tuple[Path, Path, Path]) -> bool:
    path, campaign_root, final_root = task
    value = json.loads(path.read_text())
    rewritten = repath_value(value, campaign_root, final_root)
    if rewritten == value:
        return False
    _write_json_atomic(path, rewritten)
    return True


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    metadata_dir = args.metadata_dir.resolve()
    campaign_root = args.campaign_root.resolve()
    final_root = args.final_root.resolve()
    paths = sorted(metadata_dir.glob("*.json"))
    if not paths:
        raise SystemExit(f"No metadata JSON files found in {metadata_dir}")
    if args.workers <= 0:
        raise SystemExit("--workers must be positive")
    tasks = [(path, campaign_root, final_root) for path in paths]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        changed = sum(executor.map(_rewrite_one, tasks, chunksize=16))
    print(f"Rewrote staged paths in {changed} of {len(paths)} metadata files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
