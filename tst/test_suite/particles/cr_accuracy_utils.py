"""Shared helpers for CR tracer accuracy regression tests."""

import math
import shutil
import sys
from pathlib import Path
from typing import Dict
from typing import Iterable
from typing import List
from typing import Sequence
from typing import Tuple


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts" / "particles"))
import cr_tracer_inspect  # noqa: E402


PARTICLE_OUTPUT_DIRS = (
    "ppd", "prst", "df", "dxh", "drh", "dparh", "pmom", "pspec", "pspec2",
    "psamp", "trk")


def clean_particle_outputs() -> None:
    """Remove particle outputs from the current test run directory."""
    for directory in PARTICLE_OUTPUT_DIRS:
        shutil.rmtree(directory, ignore_errors=True)
    for pattern in ("*.vtk", "*.pvtk", "*.xdmf"):
        for file_path in Path(".").glob(pattern):
            file_path.unlink()


def latest_restart_summary(path: Path = Path("prst")) -> Dict:
    """Read the latest particle restart summary."""
    return cr_tracer_inspect.summarize_restart(path)


def single_particle(summary: Dict) -> Dict:
    """Return the only particle in a one-particle restart summary."""
    particles = [
        particle
        for restart in summary["restart"]
        for particle in restart["particles"]
    ]
    if len(particles) != 1:
        raise ValueError(f"Expected one particle, found {len(particles)}")
    return particles[0]


def particles_by_species_tag(summary: Dict) -> Dict[Tuple[int, int], Dict]:
    """Return a species/tag keyed map for all restart particles."""
    result = {}
    for restart in summary["restart"]:
        for particle in restart["particles"]:
            key = (particle["species"], particle["tag"])
            if key in result:
                raise ValueError(f"Duplicate particle key {key}")
            result[key] = particle
    return result


def particle_position(particle: Dict) -> Tuple[float, float, float]:
    """Return particle position from an inspector particle record."""
    return tuple(particle["reals"][0:3])


def particle_velocity(particle: Dict) -> Tuple[float, float, float]:
    """Return particle velocity from an inspector particle record."""
    return tuple(particle["reals"][3:6])


def particle_bfield(particle: Dict) -> Tuple[float, float, float]:
    """Return sampled particle magnetic field from an inspector particle record."""
    return tuple(particle["reals"][7:10])


def norm3(values: Sequence[float]) -> float:
    """Euclidean norm of a three-vector."""
    return math.sqrt(sum(value*value for value in values))


def vector_error(a: Sequence[float], b: Sequence[float]) -> float:
    """Euclidean norm of a-b for three-vectors."""
    return norm3([a[i] - b[i] for i in range(3)])


def speed(values: Sequence[float]) -> float:
    """Magnitude of a velocity vector."""
    return norm3(values)


def gyro_exact_velocity(time: float, vperp: float = 1.0,
                        vpar: float = 0.0) -> Tuple[float, float, float]:
    """Exact velocity for the branch's uniform Bz=1 Boris sign convention."""
    return (vperp*math.cos(time), -vperp*math.sin(time), vpar)


def gyro_exact_position(time: float, vperp: float = 1.0,
                        vpar: float = 0.0) -> Tuple[float, float, float]:
    """Exact unwrapped position for x0=0 and v=(vperp,0,vpar), Bz=1."""
    return (vperp*math.sin(time), vperp*(math.cos(time) - 1.0), vpar*time)


def phase_error(vx: float, vy: float, time: float) -> float:
    """Wrapped gyro phase error for Bz=1 and initial velocity along +x."""
    angle = math.atan2(vy, vx)
    target = -time
    return abs(math.atan2(math.sin(angle - target), math.cos(angle - target)))


def fit_log_slope(x_values: Iterable[float], y_values: Iterable[float]) -> float:
    """Least-squares slope of log(y) versus log(x)."""
    xs = [math.log(float(x)) for x in x_values]
    ys = [math.log(float(y)) for y in y_values]
    n = len(xs)
    if n != len(ys) or n < 2:
        raise ValueError("fit_log_slope needs matching arrays of length >= 2")
    xbar = sum(xs)/n
    ybar = sum(ys)/n
    denom = sum((x - xbar)**2 for x in xs)
    if denom == 0.0:
        raise ValueError("fit_log_slope received identical x values")
    return sum((xs[i] - xbar)*(ys[i] - ybar) for i in range(n))/denom


def exact_bfield(profile: str, x: float, y: float, z: float,
                 b0: Sequence[float], bgrad: float = 1.0,
                 bamp: float = 0.05, bwave: float = 1.0) -> Tuple[float, float, float]:
    """Analytic B-field formulas matching part_random accuracy profiles."""
    bx, by, bz = b0
    if profile == "uniform":
        return bx, by, bz
    if profile == "linear_cross":
        return bx + bgrad*y, by + bgrad*x, bz
    if profile == "mirror":
        return bx - bgrad*x*z, by - bgrad*y*z, bz + bgrad*z*z
    if profile == "gradb":
        return bx, by, bz + bgrad*x
    if profile in ("sinusoidal_divb_free", "turbulent"):
        two_pi = 2.0*math.pi
        k = two_pi*bwave
        bx += bamp*k*math.sin(k*x)*(math.cos(k*y) - math.cos(k*z))
        by += bamp*k*math.sin(k*y)*(math.cos(k*z) - math.cos(k*x))
        bz += bamp*k*math.sin(k*z)*(math.cos(k*x) - math.cos(k*y))
        if profile == "turbulent":
            k2 = 2.0*k
            amp2 = 0.35*bamp
            bx += amp2*k2*math.sin(k2*(x + 0.13))*(
                math.cos(k2*(y - 0.07)) - math.cos(k2*(z + 0.11)))
            by += amp2*k2*math.sin(k2*(y - 0.07))*(
                math.cos(k2*(z + 0.11)) - math.cos(k2*(x + 0.13)))
            bz += amp2*k2*math.sin(k2*(z + 0.11))*(
                math.cos(k2*(x + 0.13)) - math.cos(k2*(y - 0.07)))
        return bx, by, bz
    raise ValueError(f"Unknown B profile {profile}")


def field_sample_error(summary: Dict, profile: str, b0: Sequence[float],
                       bgrad: float = 1.0, bamp: float = 0.05,
                       bwave: float = 1.0) -> float:
    """Return RMS sampled-field error for particles in a restart summary."""
    errors: List[float] = []
    for restart in summary["restart"]:
        for particle in restart["particles"]:
            x, y, z = particle_position(particle)
            exact = exact_bfield(profile, x, y, z, b0, bgrad, bamp, bwave)
            errors.append(vector_error(particle_bfield(particle), exact))
    if not errors:
        raise ValueError("No particles found for field_sample_error")
    return math.sqrt(sum(error*error for error in errors)/len(errors))
