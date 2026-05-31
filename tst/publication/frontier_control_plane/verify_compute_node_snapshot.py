#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Verify immutable submission dependencies before a compute-node run."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
from pathlib import Path

from control_plane_common import read_json, require_not_symlink, require_read_only, sha256
from control_plane_common import verify_snapshot_files


def verify(
    manifest_path: Path,
    *,
    submission_id: str | None = None,
    reservation_id: str | None = None,
    manifest_sha256: str | None = None,
) -> None:
    if manifest_sha256 is not None and sha256(manifest_path) != manifest_sha256:
        raise ValueError("Pre-submit manifest checksum differs from scheduled digest")
    require_read_only(manifest_path)
    manifest = read_json(manifest_path)
    if submission_id is not None and manifest.get("submission_id") != submission_id:
        raise ValueError("Submission ID does not match pre-submit manifest")
    if reservation_id is not None:
        reservation_path = manifest_path.parent / "reservation_id.txt"
        require_not_symlink(reservation_path)
        if not reservation_path.is_file():
            raise ValueError("Missing immutable reservation attachment")
        if reservation_path.read_bytes() != (reservation_id + "\n").encode("utf-8"):
            raise ValueError("Reservation ID does not match attachment")
        require_read_only(reservation_path)
    if manifest_sha256 is not None:
        manifest_digest_path = manifest_path.parent / "manifest_sha256.txt"
        require_not_symlink(manifest_digest_path)
        if not manifest_digest_path.is_file():
            raise ValueError("Missing immutable manifest checksum attachment")
        if manifest_digest_path.read_bytes() != (manifest_sha256 + "\n").encode("utf-8"):
            raise ValueError("Manifest checksum does not match immutable attachment")
        require_read_only(manifest_digest_path)
    verify_snapshot_files(manifest, root=Path(str(manifest["pic_root"])))
    for record in manifest["snapshot_files"]:
        require_read_only(Path(str(record["path"])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--submission-id")
    parser.add_argument("--reservation-id")
    parser.add_argument("--manifest-sha256")
    args = parser.parse_args()
    verify(
        args.manifest,
        submission_id=args.submission_id,
        reservation_id=args.reservation_id,
        manifest_sha256=args.manifest_sha256,
    )


if __name__ == "__main__":
    main()
