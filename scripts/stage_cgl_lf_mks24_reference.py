#!/usr/bin/env python3
"""Stage the pinned MKS24 arXiv source bundle with checksum provenance."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
import urllib.request


ROOT_DIR = Path(__file__).resolve().parents[1]
ARXIV_ID = "2405.02418v2"
TITLE = "Self-organization in collisionless, high-beta turbulence"
RECORD_URL = f"https://arxiv.org/abs/{ARXIV_ID}"
SOURCE_URL = f"https://arxiv.org/src/{ARXIV_ID}"
LICENSE_URL = "https://arxiv.org/licenses/nonexclusive-distrib/1.0/"
DEFAULT_OUTPUT = (
    ROOT_DIR / "build-cgl-implementation" / "cgl_lf_reference" / f"arXiv-{ARXIV_ID}"
)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_source(destination: Path) -> str:
    """Download the pinned arXiv source archive and return the final URL."""

    request = urllib.request.Request(
        SOURCE_URL,
        headers={"User-Agent": "AthenaK-CGL-LF-reference-staging/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        with destination.open("wb") as stream:
            shutil.copyfileobj(response, stream)
        return response.geturl()


def extract_archive(archive: Path, destination: Path) -> list[Path]:
    """Extract regular archive entries after rejecting unsafe paths or links."""

    destination.mkdir(parents=True)
    root = destination.resolve()
    with tarfile.open(archive, "r:*") as source:
        members = source.getmembers()
        for member in members:
            path = (destination / member.name).resolve()
            try:
                path.relative_to(root)
            except ValueError as error:
                raise ValueError(
                    f"archive path escapes destination: {member.name}"
                ) from error
            if member.issym() or member.islnk():
                raise ValueError(f"archive link is not accepted: {member.name}")
            if not member.isfile() and not member.isdir():
                raise ValueError(f"unsupported archive entry: {member.name}")
        source.extractall(destination, members=members)
    return sorted(path for path in destination.rglob("*") if path.is_file())


def parser() -> argparse.ArgumentParser:
    """Build command-line argument parsing."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Immutable staging directory to create.",
    )
    command.add_argument(
        "--archive",
        type=Path,
        help="Use an already downloaded arXiv source archive instead of downloading.",
    )
    command.add_argument(
        "--expected-sha256",
        help="Fail unless the source archive has this SHA-256 digest.",
    )
    command.add_argument(
        "--replace",
        action="store_true",
        help="Replace an existing staging directory.",
    )
    return command


def main() -> int:
    """Stage a source archive and write its machine-readable provenance."""

    args = parser().parse_args()
    target = args.output_dir.resolve()
    if target.exists() and not args.replace:
        raise ValueError(f"staging directory already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=target.parent, prefix=".mks24-stage-") as temp:
        stage = Path(temp) / target.name
        stage.mkdir()
        archive = stage / f"arXiv-{ARXIV_ID}-source.tar.gz"
        if args.archive is None:
            resolved_source_url = download_source(archive)
            source_method = "download"
        else:
            supplied_archive = args.archive.resolve()
            if not supplied_archive.is_file():
                raise ValueError(f"source archive does not exist: {supplied_archive}")
            shutil.copy2(supplied_archive, archive)
            resolved_source_url = str(supplied_archive)
            source_method = "provided_archive"

        archive_digest = sha256(archive)
        if args.expected_sha256 and archive_digest != args.expected_sha256.lower():
            raise ValueError(
                f"source archive SHA-256 mismatch: got {archive_digest}, "
                f"expected {args.expected_sha256.lower()}"
            )

        source_dir = stage / "source"
        extracted = extract_archive(archive, source_dir)
        tex_matches = sorted(source_dir.rglob("MKS24.tex"))
        if len(tex_matches) != 1:
            raise ValueError(
                f"expected one extracted MKS24.tex file; found {len(tex_matches)}"
            )

        files = {
            str(path.relative_to(stage)): sha256(path)
            for path in extracted
        }
        manifest = {
            "schema_version": 1,
            "staged_utc": datetime.now(timezone.utc).isoformat(),
            "paper": {
                "arxiv_id": ARXIV_ID,
                "title": TITLE,
                "record_url": RECORD_URL,
                "source_url": SOURCE_URL,
                "license_url": LICENSE_URL,
            },
            "source_method": source_method,
            "resolved_source": resolved_source_url,
            "archive": {
                "path": archive.name,
                "sha256": archive_digest,
            },
            "source_tex": str(tex_matches[0].relative_to(stage)),
            "extracted_file_sha256": files,
        }
        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if target.exists():
            shutil.rmtree(target)
        stage.rename(target)

    print(f"Staged MKS24 reference bundle: {target}")
    print(f"Archive SHA-256: {archive_digest}")
    print(f"Source TeX: {target / manifest['source_tex']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, tarfile.TarError) as error:
        print(f"MKS24 reference staging failed: {error}", file=sys.stderr)
        raise SystemExit(1)
