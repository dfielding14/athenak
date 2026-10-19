#!/usr/bin/env python3
"""Encode descriptive rebuilt-PDF PNG sequences as QuickTime-compatible MP4s."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

from PIL import Image


TIME_TOKEN_RE = re.compile(r"^t\d{8}\.\d{3}Myr$")
PNG_NAME_RE = re.compile(
    r"^(?P<view>.+)\."
    r"(?P<simulation>res_(?:4|8)pc(?:_(?:highmdot|lowmdot|uniform))?)\."
    r"(?P<time>t\d{8}\.\d{3}Myr)$"
)


def _write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def discover_groups(png_dir: Path) -> Dict[Tuple[str, str], List[Path]]:
    groups: Dict[Tuple[str, str], List[Path]] = defaultdict(list)
    for path in sorted(png_dir.glob("*.png")):
        match = PNG_NAME_RE.match(path.stem)
        if match is None:
            raise ValueError(f"PNG lacks semantic view/simulation/time name: {path}")
        view = match.group("view")
        simulation = match.group("simulation")
        time_token = match.group("time")
        if not TIME_TOKEN_RE.match(time_token):
            raise ValueError(f"PNG lacks movie-ready time token: {path}")
        if re.search(r"(?:^|[._-])output\d+(?:$|[._-])", path.name):
            raise ValueError(f"PNG name contains forbidden outputXX label: {path}")
        groups[(view, simulation)].append(path)

    for key, paths in groups.items():
        paths.sort(key=lambda path: PNG_NAME_RE.match(path.stem).group("time"))
        tokens = [PNG_NAME_RE.match(path.stem).group("time") for path in paths]
        if len(tokens) != len(set(tokens)):
            raise ValueError(f"Duplicate time token in movie group {key}: {tokens}")
    return dict(groups)


def _even_ceiling(value: float) -> int:
    return int(2 * math.ceil(value / 2.0))


def _movie_canvas(paths: Sequence[Path], max_dimension: int) -> Tuple[int, int]:
    dimensions = []
    for path in paths:
        with Image.open(path) as image:
            dimensions.append(image.size)
    maximum_width = max(width for width, _ in dimensions)
    maximum_height = max(height for _, height in dimensions)
    scale = min(1.0, max_dimension / max(maximum_width, maximum_height))
    return (
        _even_ceiling(maximum_width * scale),
        _even_ceiling(maximum_height * scale),
    )


def _probe_movie(ffprobe: Path, path: Path) -> Dict[str, Any]:
    command = [
        str(ffprobe),
        "-v",
        "error",
        "-select_streams",
        "v:0",
        "-show_entries",
        "stream=codec_name,profile,pix_fmt,width,height,r_frame_rate,nb_frames",
        "-show_entries",
        "format=format_name,duration,size",
        "-of",
        "json",
        str(path),
    ]
    value = json.loads(subprocess.check_output(command, text=True))
    stream = value["streams"][0]
    if stream.get("codec_name") != "h264":
        raise ValueError(f"{path} is not H.264: {stream}")
    if stream.get("pix_fmt") != "yuv420p":
        raise ValueError(f"{path} is not YUV420p: {stream}")
    if "mp4" not in value["format"].get("format_name", ""):
        raise ValueError(f"{path} is not MP4: {value['format']}")
    return value


def encode_group(task: Mapping[str, Any]) -> Dict[str, Any]:
    paths = [Path(value) for value in task["paths"]]
    output = Path(task["output"])
    ffmpeg = Path(task["ffmpeg"])
    ffprobe = Path(task["ffprobe"])
    fps = int(task["fps"])
    minimum_frames = max(len(paths), int(math.ceil(float(task["minimum_duration"]) * fps)))
    canvas_width, canvas_height = _movie_canvas(paths, int(task["max_dimension"]))
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not task["overwrite"]:
        raise FileExistsError(f"Movie already exists: {output}")

    with tempfile.TemporaryDirectory(prefix="gotham_pdf_movie_") as temporary_name:
        temporary = Path(temporary_name)
        sequence = list(paths)
        sequence.extend([paths[-1]] * (minimum_frames - len(paths)))
        for index, path in enumerate(sequence):
            (temporary / f"frame{index:06d}.png").symlink_to(path)

        filter_graph = (
            f"scale=w={canvas_width}:h={canvas_height}:"
            "force_original_aspect_ratio=decrease:flags=lanczos,"
            f"pad={canvas_width}:{canvas_height}:(ow-iw)/2:(oh-ih)/2:color=white,"
            "setsar=1,format=yuv420p"
        )
        partial = output.with_name(output.name + ".partial.mp4")
        partial.unlink(missing_ok=True)
        command = [
            str(ffmpeg),
            "-hide_banner",
            "-loglevel",
            "error",
            "-framerate",
            str(fps),
            "-i",
            str(temporary / "frame%06d.png"),
            "-vf",
            filter_graph,
            "-an",
            "-c:v",
            "libx264",
            "-preset",
            str(task["preset"]),
            "-crf",
            str(task["crf"]),
            "-profile:v",
            "high",
            "-level:v",
            "5.0",
            "-pix_fmt",
            "yuv420p",
            "-movflags",
            "+faststart",
            "-threads",
            str(task["ffmpeg_threads"]),
            "-y",
            str(partial),
        ]
        subprocess.run(command, check=True)
        partial.replace(output)

    probe = _probe_movie(ffprobe, output)
    return {
        "view": task["view"],
        "simulation": task["simulation"],
        "output": str(output),
        "source_frame_count": len(paths),
        "encoded_frame_count": minimum_frames,
        "first_frame": str(paths[0]),
        "last_frame": str(paths[-1]),
        "canvas": [canvas_width, canvas_height],
        "probe": probe,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("png_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ffmpeg", type=Path, required=True)
    parser.add_argument("--ffprobe", type=Path, required=True)
    parser.add_argument("--fps", type=int, default=5)
    parser.add_argument("--minimum-duration", type=float, default=3.0)
    parser.add_argument("--max-dimension", type=int, default=2160)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--ffmpeg-threads", type=int, default=4)
    parser.add_argument("--preset", default="slow")
    parser.add_argument("--crf", type=int, default=20)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    png_dir = args.png_dir.resolve()
    output_dir = args.output_dir.resolve()
    ffmpeg = args.ffmpeg.resolve()
    ffprobe = args.ffprobe.resolve()
    if not png_dir.is_dir():
        raise SystemExit(f"PNG directory does not exist: {png_dir}")
    for executable in (ffmpeg, ffprobe):
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise SystemExit(f"Required executable is unavailable: {executable}")
    if args.fps <= 0 or args.workers <= 0 or args.ffmpeg_threads <= 0:
        raise SystemExit("FPS and worker counts must be positive")
    if args.minimum_duration <= 0.0 or args.max_dimension <= 0:
        raise SystemExit("Minimum duration and maximum dimension must be positive")

    groups = discover_groups(png_dir)
    tasks = []
    for (view, simulation), paths in sorted(groups.items()):
        output = output_dir / f"{view}.{simulation}.mp4"
        tasks.append(
            {
                "view": view,
                "simulation": simulation,
                "paths": [str(path) for path in paths],
                "output": str(output),
                "ffmpeg": str(ffmpeg),
                "ffprobe": str(ffprobe),
                "fps": args.fps,
                "minimum_duration": args.minimum_duration,
                "max_dimension": args.max_dimension,
                "ffmpeg_threads": args.ffmpeg_threads,
                "preset": args.preset,
                "crf": args.crf,
                "overwrite": args.overwrite,
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    results: List[Dict[str, Any]] = []
    failures: List[Dict[str, str]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_tasks = {executor.submit(encode_group, task): task for task in tasks}
        for future in concurrent.futures.as_completed(future_tasks):
            task = future_tasks[future]
            try:
                result = future.result()
            except Exception as exc:
                failures.append(
                    {
                        "view": str(task["view"]),
                        "simulation": str(task["simulation"]),
                        "error": str(exc),
                    }
                )
                print(f"FAIL {task['view']}.{task['simulation']}: {exc}")
                continue
            results.append(result)
            print(
                f"DONE {result['view']}.{result['simulation']}: "
                f"{result['source_frame_count']} source frames"
            )

    results.sort(key=lambda value: (value["view"], value["simulation"]))
    manifest = {
        "schema_version": 1,
        "kind": "rebuilt_pdf_movie_campaign",
        "png_dir": str(png_dir),
        "output_dir": str(output_dir),
        "ffmpeg": str(ffmpeg),
        "ffprobe": str(ffprobe),
        "fps": args.fps,
        "minimum_duration": args.minimum_duration,
        "max_dimension": args.max_dimension,
        "preset": args.preset,
        "crf": args.crf,
        "movie_count": len(results),
        "failed_movie_count": len(failures),
        "source_frame_count": sum(item["source_frame_count"] for item in results),
        "movies": results,
        "failures": failures,
    }
    _write_json_atomic(output_dir / "movie_manifest.json", manifest)
    print(
        f"Wrote {len(results)} movies from {manifest['source_frame_count']} source "
        f"frames to {output_dir}"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
