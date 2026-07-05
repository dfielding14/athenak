#!/usr/bin/env python3
"""Analyze the CR pusher accuracy-test sweep.

The companion runner writes one JSON record per case and one tracked-particle
file per run.  This analyzer turns those tracks into the table requested by
docs/source/modules/cr_pusher_accuracy_test_plan.md.
"""

import argparse
import csv
import json
import math
import re
import struct
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


Vec3 = Tuple[float, float, float]
ParticleKey = Tuple[int, int]
TrackRecord = Tuple[float, ...]
RICH_FIELDS = ("tag", "time", "x", "y", "z", "vx", "vy", "vz",
               "bx", "by", "bz", "k1", "k2", "k3", "db1", "db2", "db3",
               "jmag")
RICH_RECORD_FIELDS = RICH_FIELDS[2:]
RICH_DIAGNOSTIC_SLICE = slice(6, 16)
TRACK_HEADER_MARKER = b"# AthenaK tracked particle data at time="


class TrackFrame:
    def __init__(self, time: float, cycle: int, ntracked: int,
                 ntrack_per_species: int, track_per_species: bool,
                 particles: Dict[ParticleKey, TrackRecord],
                 nfields: int = 6,
                 fields: Tuple[str, ...] = ()):
        self.time = time
        self.cycle = cycle
        self.ntracked = ntracked
        self.ntrack_per_species = ntrack_per_species
        self.track_per_species = track_per_species
        self.particles = particles
        self.nfields = nfields
        self.fields = fields


def _norm(values: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return sum(a[i] * b[i] for i in range(3))


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _vadd(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _vsub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vmul(a: Sequence[float], scalar: float) -> Vec3:
    return (a[0] * scalar, a[1] * scalar, a[2] * scalar)


def _periodic_delta(actual: Sequence[float], expected: Sequence[float],
                    box: Sequence[float]) -> Vec3:
    delta = [actual[i] - expected[i] for i in range(3)]
    for i, length in enumerate(box):
        if length > 0.0:
            delta[i] -= round(delta[i] / length) * length
    return (delta[0], delta[1], delta[2])


def _position(record: Sequence[float]) -> Vec3:
    return (record[0], record[1], record[2])


def _velocity(record: Sequence[float]) -> Vec3:
    return (record[3], record[4], record[5])


def _finite_record(record: Sequence[float]) -> bool:
    return all(math.isfinite(value) for value in record)


def _mean_square(values: Iterable[float]) -> float:
    values = list(values)
    if not values:
        return math.nan
    return sum(value * value for value in values) / len(values)


def _rms(values: Iterable[float]) -> float:
    return math.sqrt(_mean_square(values))


def _mass(species: int, meta: Dict) -> float:
    return float(meta["min_mass"]) * float(meta["mass_log_spacing"]) ** species


def _box_lengths(meta: Dict) -> Vec3:
    return (
        float(meta.get("x1max", 0.5)) - float(meta.get("x1min", -0.5)),
        float(meta.get("x2max", 0.5)) - float(meta.get("x2min", -0.5)),
        float(meta.get("x3max", 0.5)) - float(meta.get("x3min", -0.5)),
    )


def _exact_bfield(meta: Dict, x: float, y: float, z: float) -> Vec3:
    bx = float(meta.get("B0x", 0.0))
    by = float(meta.get("B0y", 0.0))
    bz = float(meta.get("B0z", 1.0))
    profile = str(meta.get("B_profile", "uniform"))
    bgrad = float(meta.get("Bgrad", 1.0))
    bamp = float(meta.get("Bamp", 0.05))
    bwave = float(meta.get("Bwave_number", 1.0))

    if profile == "uniform":
        return bx, by, bz
    if profile == "linear_cross":
        return bx + bgrad * y, by + bgrad * x, bz
    if profile == "mirror":
        return bx - bgrad * x * z, by - bgrad * y * z, bz + bgrad * z * z
    if profile == "gradb":
        return bx, by, bz + bgrad * x
    if profile in ("sinusoidal_divb_free", "turbulent"):
        k = 2.0 * math.pi * bwave
        bx += bamp * k * math.sin(k * x) * (math.cos(k * y) - math.cos(k * z))
        by += bamp * k * math.sin(k * y) * (math.cos(k * z) - math.cos(k * x))
        bz += bamp * k * math.sin(k * z) * (math.cos(k * x) - math.cos(k * y))
        if profile == "turbulent":
            k2 = 2.0 * k
            amp2 = 0.35 * bamp
            bx += amp2 * k2 * math.sin(k2 * (x + 0.13)) * (
                math.cos(k2 * (y - 0.07)) - math.cos(k2 * (z + 0.11)))
            by += amp2 * k2 * math.sin(k2 * (y - 0.07)) * (
                math.cos(k2 * (z + 0.11)) - math.cos(k2 * (x + 0.13)))
            bz += amp2 * k2 * math.sin(k2 * (z + 0.11)) * (
                math.cos(k2 * (x + 0.13)) - math.cos(k2 * (y - 0.07)))
        return bx, by, bz
    raise ValueError(f"Unknown B_profile={profile}")


def _magnetic_moment(velocity: Vec3, bfield: Vec3) -> float:
    bmag = _norm(bfield)
    if bmag == 0.0:
        return math.nan
    speed2 = _dot(velocity, velocity)
    vpar = _dot(velocity, bfield) / bmag
    return (speed2 - vpar * vpar) / bmag


def _exact_uniform_state(pos0: Vec3, vel0: Vec3, mass: float, bfield: Vec3,
                         dt: float) -> Tuple[Vec3, Vec3]:
    bmag = _norm(bfield)
    if bmag == 0.0 or mass <= 0.0:
        return _vadd(pos0, _vmul(vel0, dt)), vel0
    bhat = _vmul(bfield, 1.0 / bmag)
    omega = bmag / mass
    angle = omega * dt
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    vpar = _vmul(bhat, _dot(vel0, bhat))
    vperp = _vsub(vel0, vpar)
    b_cross_v = _cross(bhat, vperp)

    vel = _vadd(vpar, _vsub(_vmul(vperp, cos_a), _vmul(b_cross_v, sin_a)))
    dpos = _vadd(
        _vmul(vpar, dt),
        _vadd(_vmul(vperp, sin_a / omega),
              _vmul(b_cross_v, (cos_a - 1.0) / omega)),
    )
    return _vadd(pos0, dpos), vel


def _phase_error(actual: Vec3, initial: Vec3, mass: float, bfield: Vec3,
                 dt: float) -> float:
    if abs(bfield[0]) > 1.0e-14 or abs(bfield[1]) > 1.0e-14:
        return math.nan
    omega = abs(bfield[2]) / mass
    initial_angle = math.atan2(initial[1], initial[0])
    actual_angle = math.atan2(actual[1], actual[0])
    target = initial_angle - math.copysign(omega * dt, bfield[2])
    return abs(math.atan2(math.sin(actual_angle - target),
                          math.cos(actual_angle - target)))


def read_trk_file(path: Path) -> List[TrackFrame]:
    blob = path.read_bytes()
    time_pattern = re.compile(
        r"# AthenaK tracked particle data at time=\s*([0-9.eE+-]+)")
    key_pattern = re.compile(r"([A-Za-z_]+)=\s*([^ \t\n]+)")
    frames: List[TrackFrame] = []
    offset = 0
    while True:
        start = blob.find(TRACK_HEADER_MARKER, offset)
        if start < 0:
            break
        cursor = start
        header_lines: List[str] = []
        while True:
            line_end = blob.find(b"\n", cursor)
            if line_end < 0:
                raise ValueError(f"{path} has an unterminated track header")
            raw_line = blob[cursor:line_end]
            cursor = line_end + 1
            if raw_line.strip() == b"":
                break
            header_lines.append(raw_line.decode("ascii", errors="replace"))

        header_text = " ".join(header_lines)
        time_match = time_pattern.search(header_text)
        if time_match is None:
            raise ValueError(f"{path} has an unrecognized track header: {header_text}")
        header_values = {key: value for key, value in key_pattern.findall(header_text)}

        time = float(time_match.group(1))
        cycle = int(header_values["cycle"])
        ntracked = int(header_values["ntracked_prtcls"])
        ntrack_per_species = int(header_values.get("ntrack_per_species", ntracked))
        track_per_species = bool(int(header_values.get("track_per_species", "1")))
        record_count = int(header_values.get("record_count", str(ntracked)))
        if "nfields" in header_values:
            nfields = int(header_values["nfields"])
        else:
            nfields = _infer_track_nfields(path, blob, cursor, record_count)

        if nfields not in (6, 18):
            raise ValueError(f"{path} has unsupported trk nfields={nfields}")
        fields = tuple(header_values.get(
            "fields", ",".join(RICH_FIELDS if nfields == 18 else ())).split(","))
        if fields == ("",):
            fields = ()
        if nfields == 18 and fields and fields != RICH_FIELDS:
            raise ValueError(f"{path} has unsupported rich trk fields={fields}")
        payload_start = cursor
        payload_bytes = nfields * record_count * struct.calcsize("<f")
        payload_end = payload_start + payload_bytes
        if payload_end > len(blob):
            raise ValueError(f"{path} has an incomplete track payload at t={time}")
        values = struct.unpack_from("<" + str(nfields * record_count) + "f", blob,
                                    payload_start)
        particles: Dict[ParticleKey, TrackRecord] = {}
        for index in range(record_count):
            base = nfields * index
            if nfields == 18:
                record_tag = float(values[base])
                record_time = float(values[base + 1])
                if not math.isfinite(record_tag) or not math.isfinite(record_time):
                    raise ValueError(f"{path} has non-finite rich tag/time at t={time}")
                output_tag = int(round(record_tag))
                record = tuple(float(value) for value in values[base + 2:base + 18])
                if len(record) != len(RICH_RECORD_FIELDS):
                    raise ValueError(f"{path} has malformed rich record at t={time}")
                if not _finite_record(record[RICH_DIAGNOSTIC_SLICE]):
                    raise ValueError(
                        f"{path} has non-finite rich diagnostics at t={time}")
            else:
                output_tag = index
                record = tuple(float(value) for value in values[base:base + 6])
            if track_per_species:
                species = output_tag // ntrack_per_species
                tag = output_tag % ntrack_per_species
            else:
                species = 0
                tag = output_tag
            particles[(species, tag)] = record
        frames.append(TrackFrame(time, cycle, ntracked, ntrack_per_species,
                                 track_per_species, particles, nfields, fields))
        offset = payload_end
    if not frames:
        raise ValueError(f"No tracked-particle frames found in {path}")
    return frames


def _track_payload_boundary_ok(blob: bytes, payload_end: int) -> bool:
    cursor = payload_end
    while cursor < len(blob) and blob[cursor:cursor + 1] in b" \t\r\n":
        cursor += 1
    return cursor == len(blob) or blob.startswith(TRACK_HEADER_MARKER, cursor)


def _infer_track_nfields(path: Path, blob: bytes, payload_start: int,
                         record_count: int) -> int:
    for nfields in (18, 6):
        payload_end = payload_start + nfields * record_count * struct.calcsize("<f")
        if payload_end <= len(blob) and _track_payload_boundary_ok(blob, payload_end):
            return nfields
    raise ValueError(f"{path} is missing nfields and payload size is not 6- or 18-field")


def _case_track_files(run_dir: Path) -> List[Path]:
    trk_dir = run_dir / "trk"
    files = sorted(trk_dir.glob("*.trk"))
    if not files:
        files = sorted(trk_dir.glob("*/*.trk"))
    if not files:
        raise FileNotFoundError(f"No .trk file under {trk_dir}")
    return files


def _read_case_tracks(run_dir: Path) -> List[TrackFrame]:
    merged: Dict[Tuple[float, int], TrackFrame] = {}
    for path in _case_track_files(run_dir):
        for frame in read_trk_file(path):
            key = (round(frame.time, 12), frame.cycle)
            if key not in merged:
                merged[key] = TrackFrame(frame.time, frame.cycle, frame.ntracked,
                                         frame.ntrack_per_species,
                                         frame.track_per_species, {},
                                         frame.nfields, frame.fields)
            target = merged[key]
            if (target.ntracked != frame.ntracked or
                    target.ntrack_per_species != frame.ntrack_per_species or
                    target.track_per_species != frame.track_per_species):
                raise ValueError(f"Inconsistent trk metadata for frame t={frame.time}")
            if target.nfields != frame.nfields:
                raise ValueError(f"Inconsistent nfields for frame t={frame.time}")
            duplicates = sorted(set(target.particles) & set(frame.particles))
            if duplicates:
                raise ValueError(
                    f"Duplicate tracked particle keys at t={frame.time}: {duplicates[:3]}")
            target.particles.update(frame.particles)
    return [merged[key] for key in sorted(merged)]


def _read_time_file(run_dir: Path) -> Optional[float]:
    time_file = run_dir / "run.time"
    if not time_file.exists():
        return None
    for line in time_file.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] == "real":
            return float(parts[1])
    return None


def _run_had_fatal(run_dir: Path) -> bool:
    for name in ("stdout.txt", "stderr.txt"):
        path = run_dir / name
        if path.exists() and "FATAL ERROR" in path.read_text(errors="replace"):
            return True
    return False


def _load_cases(manifest: Path) -> List[Dict]:
    cases = []
    for line in manifest.read_text().splitlines():
        stripped = line.strip()
        if stripped:
            cases.append(json.loads(stripped))
    if not cases:
        raise ValueError(f"{manifest} contains no cases")
    return cases


def _frame_times(frames: Sequence[TrackFrame]) -> Dict[float, TrackFrame]:
    return {round(frame.time, 12): frame for frame in frames}


def _species_list(meta: Dict) -> List[int]:
    return list(range(int(meta["nspecies"])))


def _empty_species_stats(meta: Dict) -> Dict[int, Dict[str, List[float]]]:
    return {
        species: {
            "dx": [],
            "dv": [],
            "phase": [],
            "dmu": [],
            "de": [],
        }
        for species in _species_list(meta)
    }


def _add_invariant_errors(stats: Dict[int, Dict[str, List[float]]],
                          frames: Sequence[TrackFrame], meta: Dict) -> None:
    first = frames[0]
    initial = first.particles
    initial_mu: Dict[ParticleKey, float] = {}
    initial_e: Dict[ParticleKey, float] = {}
    for key, record in initial.items():
        if not _finite_record(record):
            continue
        pos = _position(record)
        vel = _velocity(record)
        bfield = _exact_bfield(meta, *pos)
        initial_mu[key] = _magnetic_moment(vel, bfield)
        initial_e[key] = 0.5 * _dot(vel, vel)

    for frame in frames:
        for key, record in frame.particles.items():
            if not _finite_record(record) or key not in initial_mu:
                continue
            species, _ = key
            pos = _position(record)
            vel = _velocity(record)
            bfield = _exact_bfield(meta, *pos)
            mu0 = initial_mu.get(key, math.nan)
            e0 = initial_e.get(key, math.nan)
            mu = _magnetic_moment(vel, bfield)
            energy = 0.5 * _dot(vel, vel)
            if math.isfinite(mu0) and abs(mu0) > 0.0:
                stats[species]["dmu"].append(abs((mu - mu0) / mu0))
            if math.isfinite(e0) and abs(e0) > 0.0:
                stats[species]["de"].append(abs((energy - e0) / e0))


def _analyze_uniform(meta: Dict, frames: Sequence[TrackFrame]) -> Dict[int, Dict[str, List[float]]]:
    stats = _empty_species_stats(meta)
    first = frames[0]
    box = _box_lengths(meta)
    bfield = (
        float(meta.get("B0x", 0.0)),
        float(meta.get("B0y", 0.0)),
        float(meta.get("B0z", 1.0)),
    )
    for frame in frames:
        dt = frame.time - first.time
        for key, record in frame.particles.items():
            species, _ = key
            initial = first.particles[key]
            if not _finite_record(record) or not _finite_record(initial):
                continue
            mass = _mass(species, meta)
            exact_pos, exact_vel = _exact_uniform_state(
                _position(initial), _velocity(initial), mass, bfield, dt)
            dx = _norm(_periodic_delta(_position(record), exact_pos, box))
            dv = _norm(_vsub(_velocity(record), exact_vel))
            phase = _phase_error(_velocity(record), _velocity(initial), mass,
                                 bfield, dt)
            stats[species]["dx"].append(dx)
            stats[species]["dv"].append(dv)
            if math.isfinite(phase):
                stats[species]["phase"].append(phase)
    _add_invariant_errors(stats, frames, meta)
    return stats


def _compare_to_reference(meta: Dict, frames: Sequence[TrackFrame],
                          ref_frames: Sequence[TrackFrame]) -> Dict[int, Dict[str, List[float]]]:
    stats = _empty_species_stats(meta)
    box = _box_lengths(meta)
    by_time = _frame_times(frames)
    ref_by_time = _frame_times(ref_frames)
    for time_key in sorted(set(by_time) & set(ref_by_time)):
        frame = by_time[time_key]
        ref = ref_by_time[time_key]
        for key in sorted(set(frame.particles) & set(ref.particles)):
            species, _ = key
            if not _finite_record(frame.particles[key]) or not _finite_record(ref.particles[key]):
                continue
            dx = _norm(_periodic_delta(
                _position(frame.particles[key]), _position(ref.particles[key]), box))
            dv = _norm(_vsub(_velocity(frame.particles[key]),
                             _velocity(ref.particles[key])))
            stats[species]["dx"].append(dx)
            stats[species]["dv"].append(dv)
    _add_invariant_errors(stats, frames, meta)
    return stats


def _max_step_displacement(frames: Sequence[TrackFrame], meta: Dict) -> float:
    box = _box_lengths(meta)
    max_step = 0.0
    for left, right in zip(frames, frames[1:]):
        for key in sorted(set(left.particles) & set(right.particles)):
            if not _finite_record(left.particles[key]) or not _finite_record(right.particles[key]):
                continue
            delta = _periodic_delta(
                _position(right.particles[key]), _position(left.particles[key]), box)
            max_step = max(max_step, _norm(delta))
    return max_step


def _complete_track(frames: Sequence[TrackFrame], meta: Dict) -> bool:
    expected = int(meta["ntrack"]) * int(meta["nspecies"])
    for frame in frames:
        if frame.ntracked != expected or len(frame.particles) != expected:
            return False
        for record in frame.particles.values():
            if not all(math.isfinite(value) for value in record):
                return False
            if frame.nfields == 18 and len(record) != len(RICH_RECORD_FIELDS):
                return False
            if frame.nfields == 18 and not _finite_record(record[RICH_DIAGNOSTIC_SLICE]):
                return False
    return True


def _summarize_stats(stats: Dict[int, Dict[str, List[float]]]) -> Dict[int, Dict[str, float]]:
    result: Dict[int, Dict[str, float]] = {}
    for species, values in stats.items():
        dx = values["dx"]
        dv = values["dv"]
        dmu = values["dmu"]
        de = values["de"]
        phase = values["phase"]
        result[species] = {
            "max_dx": max(dx) if dx else math.nan,
            "rms_dx": _rms(dx) if dx else math.nan,
            "max_dv": max(dv) if dv else math.nan,
            "rms_dv": _rms(dv) if dv else math.nan,
            "max_phase": max(phase) if phase else math.nan,
            "dmu_mu": max(dmu) if dmu else math.nan,
            "dE_E": max(de) if de else math.nan,
        }
    return result


def _case_mode(meta: Dict) -> str:
    if meta["mode"] == "global":
        return "global"
    if meta["mode"] == "reference":
        return "reference"
    return "per_particle"


def _format_float(value: float) -> str:
    if value is None or not math.isfinite(value):
        return "nan"
    return f"{value:.8e}"


def analyze_cases(cases: Sequence[Dict]) -> Tuple[List[Dict], Dict]:
    frames_by_case: Dict[str, List[TrackFrame]] = {}
    hard_errors: List[str] = []
    for meta in cases:
        run_dir = Path(meta["run_dir"])
        try:
            frames_by_case[meta["case"]] = _read_case_tracks(run_dir)
        except Exception as exc:  # noqa: BLE001
            hard_errors.append(f"{meta['case']}: {exc}")
            frames_by_case[meta["case"]] = []
        if _run_had_fatal(run_dir):
            hard_errors.append(f"{meta['case']}: FATAL ERROR appeared in logs")

    ref_frames = {
        meta["test"]: frames_by_case[meta["case"]]
        for meta in cases
        if meta["mode"] == "reference" and frames_by_case.get(meta["case"])
    }

    rows: List[Dict] = []
    for meta in cases:
        frames = frames_by_case.get(meta["case"], [])
        if not frames:
            continue
        if not _complete_track(frames, meta):
            hard_errors.append(f"{meta['case']}: incomplete/non-finite tracked frames")
        if meta["test"] == "uniform_b":
            stats = _analyze_uniform(meta, frames)
        elif meta["test"] in ref_frames and meta["mode"] != "reference":
            stats = _compare_to_reference(meta, frames, ref_frames[meta["test"]])
        elif meta["test"] == "boundary_crossing" and "boundary_reference" in ref_frames:
            stats = _compare_to_reference(meta, frames, ref_frames["boundary_reference"])
        else:
            stats = _empty_species_stats(meta)
            _add_invariant_errors(stats, frames, meta)
        summary = _summarize_stats(stats)
        max_step = _max_step_displacement(frames, meta)
        for species, metrics in summary.items():
            rows.append({
                "test": meta["test"],
                "case": meta["case"],
                "mode": _case_mode(meta),
                "gyro_fraction": float(meta["gyro_fraction"]),
                "species": species,
                "mass": _mass(species, meta),
                "frames": len(frames),
                "runtime_s": _read_time_file(Path(meta["run_dir"])),
                "max_step_dx": max_step,
                **metrics,
            })

    checks = _acceptance_checks(cases, rows, frames_by_case, hard_errors)
    return rows, checks


def _row_by_case_species(rows: Sequence[Dict]) -> Dict[Tuple[str, int], Dict]:
    return {(row["case"], int(row["species"])): row for row in rows}


def _monotonic_non_decreasing(values: Sequence[float], rtol: float = 0.05) -> bool:
    clean = [value for value in values if math.isfinite(value)]
    if len(clean) < 2:
        return False
    for left, right in zip(clean, clean[1:]):
        if right + rtol * max(abs(left), 1.0e-300) < left:
            return False
    return True


def _acceptance_checks(cases: Sequence[Dict], rows: Sequence[Dict],
                       frames_by_case: Dict[str, List[TrackFrame]],
                       hard_errors: Sequence[str]) -> Dict:
    by_case = {meta["case"]: meta for meta in cases}
    by_case_species = _row_by_case_species(rows)
    checks: Dict[str, object] = {
        "hard_errors": list(hard_errors),
        "passed": len(hard_errors) == 0,
        "notes": [],
    }

    def add_note(text: str) -> None:
        checks["notes"].append(text)  # type: ignore[index]

    # Energy conservation should be at roundoff-ish levels for the magnetic-only
    # Boris pusher. Keep this as a hard correctness check, but loose enough for
    # float trk output.
    max_energy = max(
        (row["dE_E"] for row in rows if math.isfinite(row["dE_E"])),
        default=math.nan)
    checks["max_dE_E"] = max_energy
    if not math.isfinite(max_energy) or max_energy > 1.0e-5:
        checks["passed"] = False
        add_note(f"Energy conservation check failed: max dE/E={max_energy:.3e}")

    # Uniform g=0.05 per-particle and global paths should be close against the
    # analytic orbit. This compares their analytic RMS errors, not bit patterns.
    uniform_global = [
        row for row in rows
        if row["test"] == "uniform_b" and row["mode"] == "global" and
        abs(row["gyro_fraction"] - 0.05) < 1.0e-12
    ]
    uniform_pp = [
        row for row in rows
        if row["test"] == "uniform_b" and row["mode"] == "per_particle" and
        abs(row["gyro_fraction"] - 0.05) < 1.0e-12
    ]
    if uniform_global and uniform_pp:
        global_rms = max(row["rms_dv"] for row in uniform_global)
        pp_rms = max(row["rms_dv"] for row in uniform_pp)
        checks["uniform_g005_global_rms_dv"] = global_rms
        checks["uniform_g005_pp_rms_dv"] = pp_rms
        if pp_rms > max(2.0 * global_rms, 2.0e-5):
            add_note(
                "Uniform g=0.05 per-particle velocity error is larger than "
                "the old global-subcycle error; inspect before relaxing production.")
    else:
        add_note("Uniform g=0.05 global/per-particle comparison was incomplete.")

    # Smooth-field convergence is advisory: it can be slightly nonmonotone for a
    # short tracked subset, so report rather than fail the whole batch.
    smooth_cases = [
        meta for meta in cases
        if meta["test"] == "smooth_divb0" and meta["mode"] == "per_particle"
    ]
    smooth_cases.sort(key=lambda meta: float(meta["gyro_fraction"]))
    smooth_max = []
    for meta in smooth_cases:
        case_rows = [row for row in rows if row["case"] == meta["case"]]
        if case_rows:
            smooth_max.append(max(row["rms_dv"] for row in case_rows))
    checks["smooth_error_increases_with_gyro_fraction"] = _monotonic_non_decreasing(
        smooth_max)
    if not checks["smooth_error_increases_with_gyro_fraction"]:
        add_note("Smooth-field RMS errors were not strictly monotone in this short run.")

    boundary_rows = [row for row in rows if row["test"] == "boundary_crossing"]
    if boundary_rows:
        max_step = max(row["max_step_dx"] for row in boundary_rows)
        checks["boundary_max_step_dx"] = max_step
        if max_step > 0.35:
            checks["passed"] = False
            add_note(f"Boundary tracked-particle jump check failed: {max_step:.3e}")
    boundary_cases = [meta for meta in cases if meta["test"] == "boundary_crossing"]
    for meta in boundary_cases:
        frames = frames_by_case.get(meta["case"], [])
        if frames and len(frames[-1].particles) != int(meta["ntrack"]) * int(meta["nspecies"]):
            checks["passed"] = False
            add_note(f"{meta['case']} ended with an incomplete tracked-particle set")

    # Reference cases are allowed to have NaN dx/dv table entries; their role is
    # to provide comparison frames.
    for meta in by_case.values():
        if meta["mode"] == "reference" and not frames_by_case.get(meta["case"]):
            checks["passed"] = False
            add_note(f"Reference case {meta['case']} produced no frames")

    return checks


def write_csv(path: Path, rows: Sequence[Dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "test", "case", "mode", "gyro_fraction", "species", "mass",
        "max_dx", "rms_dx", "max_dv", "rms_dv", "max_phase", "dmu_mu",
        "dE_E", "max_step_dx", "frames", "runtime_s",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def write_markdown(path: Path, rows: Sequence[Dict], checks: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# CR Pusher Accuracy Sweep Results",
        "",
        f"Overall hard-check status: {'PASS' if checks['passed'] else 'FAIL'}",
        "",
        "## Acceptance Checks",
        "",
    ]
    for key, value in checks.items():
        if key == "notes":
            continue
        if isinstance(value, float):
            lines.append(f"- {key}: {_format_float(value)}")
        else:
            lines.append(f"- {key}: {value}")
    notes = checks.get("notes", [])
    if notes:
        lines.append("")
        lines.append("## Notes")
        lines.append("")
        for note in notes:
            lines.append(f"- {note}")

    lines.extend([
        "",
        "## Error Table",
        "",
        "| test | mode | gyro | species | max_dx | rms_dx | max_dv | rms_dv | dmu/mu | dE/E |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in rows:
        lines.append(
            "| {test} | {mode} | {gyro:.3g} | {species} | {max_dx} | {rms_dx} | "
            "{max_dv} | {rms_dv} | {dmu} | {de} |".format(
                test=row["test"],
                mode=row["mode"],
                gyro=row["gyro_fraction"],
                species=row["species"],
                max_dx=_format_float(row["max_dx"]),
                rms_dx=_format_float(row["rms_dx"]),
                max_dv=_format_float(row["max_dv"]),
                rms_dv=_format_float(row["rms_dv"]),
                dmu=_format_float(row["dmu_mu"]),
                de=_format_float(row["dE_E"]),
            ))
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--write-json", type=Path, required=True)
    parser.add_argument("--write-csv", type=Path, required=True)
    parser.add_argument("--write-md", type=Path, required=True)
    parser.add_argument("--fail-on-error", action="store_true")
    args = parser.parse_args()

    cases = _load_cases(args.manifest)
    rows, checks = analyze_cases(cases)
    write_csv(args.write_csv, rows)
    args.write_json.parent.mkdir(parents=True, exist_ok=True)
    args.write_json.write_text(json.dumps({
        "cases": cases,
        "rows": rows,
        "checks": checks,
    }, indent=2, sort_keys=True) + "\n")
    write_markdown(args.write_md, rows, checks)

    print(f"Wrote {args.write_csv}")
    print(f"Wrote {args.write_json}")
    print(f"Wrote {args.write_md}")
    if args.fail_on_error and not checks["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
