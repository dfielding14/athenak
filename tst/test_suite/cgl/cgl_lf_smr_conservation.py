"""Conservation at shared LF faces with independent transverse AMR ghost donors."""
from pathlib import Path
import struct
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis/python"))
import athena_read  # noqa: E402

INPUT_NAME = "cgl_lf_smr_conservation.athinput"
INPUT = str(ROOT / "inputs/tests" / INPUT_NAME)
BASENAME = "cgl_lf_smr_conservation"
CASES = [(2, primitive, integrator) for primitive in (False, True)
         for integrator in ("sts", "explicit")]
CASES += [(3, primitive, "sts") for primitive in (False, True)]


def flags(directory, dimension, primitive, integrator):
    result = ["-d", str(directory),
              "mesh_refinement/prolong_primitives=" + str(primitive).lower(),
              "mhd/cgl_heat_flux_integrator=" + integrator]
    if integrator == "explicit":
        result.append("time/sts_integrator=none")
    if dimension == 3:
        result += [f"mesh/nx{axis}=16" for axis in (1, 2, 3)]
        result += [f"meshblock/nx{axis}=4" for axis in (1, 2, 3)]
    return result


def read_state(path):
    """Read this fixture's double-precision MHD+one-scalar restart payload.

    Validate the serialized per-block byte count so a layout change fails loudly.
    The fixture has no other physics modules or auxiliary restart state.
    """
    raw = Path(path).read_bytes()
    header = raw.index(b"<par_end>\n") + len(b"<par_end>\n")
    nmb, root_level = struct.unpack_from("=ii", raw, header)
    root_indices = struct.unpack_from("=19i", raw, header + 8 + 72)
    indices = struct.unpack_from("=19i", raw, header + 8 + 72 + 76)
    ng, nx, ny, nz = indices[:4]
    time, dt, cycle = struct.unpack_from("=ddi", raw, header + 8 + 72 + 152)
    logical = np.frombuffer(raw, dtype="=i4", count=nmb*4,
                           offset=header + 252).reshape(nmb, 4)
    n1 = nx + 2*ng
    n2 = ny + 2*ng if ny > 1 else 1
    n3 = nz + 2*ng if nz > 1 else 1
    shapes = [(7, n3, n2, n1), (n3, n2, n1 + 1),
              (n3, n2 + 1, n1), (n3 + 1, n2, n1)]
    sizes = [int(np.prod(shape)) for shape in shapes]
    offset = len(raw) - nmb*sum(sizes)*8
    assert struct.unpack_from("=Q", raw, offset - 8)[0] == sum(sizes)*8, (
        "unsupported restart layout: this regression requires double precision MHD")
    blocks = np.frombuffer(raw, dtype="=f8", count=nmb*sum(sizes),
                           offset=offset).reshape(nmb, sum(sizes))
    arrays = []
    start = 0
    for size, shape in zip(sizes, shapes):
        arrays.append(blocks[:, start:start + size].reshape((nmb,) + shape))
        start += size
    u, b1, b2, b3 = arrays
    bcc = np.stack((.5*(b1[:, :, :, :-1] + b1[:, :, :, 1:]),
                    .5*(b2[:, :, :-1, :] + b2[:, :, 1:, :]),
                    .5*(b3[:, :-1, :, :] + b3[:, 1:, :, :])), axis=1)
    active = tuple(slice(ng, ng + n) if n > 1 else slice(0, 1)
                   for n in (nz, ny, nx))
    u = u[(slice(None), slice(None)) + active]
    bcc = bcc[(slice(None), slice(None)) + active]
    rho = u[:, 0]
    bsqr = np.sum(bcc*bcc, axis=1)
    bmag = np.sqrt(bsqr)
    internal = u[:, 4] - .5*np.sum(u[:, 1:4]**2, axis=1)/rho - .5*bsqr
    ratio = np.exp(u[:, 5]/rho - 2*np.log(rho) + 3*np.log(bmag))
    pperp = internal*ratio/(.5 + ratio)
    dimension = 3 if nz > 1 else 2
    volume = (2.0**(-dimension*(logical[:, 3] - root_level)) /
              np.prod(root_indices[1:4]))
    totals = {"energy": float(np.sum(u[:, 4].reshape(nmb, -1).sum(axis=1)*volume)),
              "mu": float(np.sum((pperp/bmag).reshape(nmb, -1).sum(axis=1)*volume))}
    assert np.all(np.isfinite(u)) and np.all(np.isfinite(pperp))
    assert np.min(rho) > 0 and np.min(pperp) > 0
    order = np.lexsort(logical.T[::-1])
    return dict(time=time, cycle=cycle, totals=totals, nmb=nmb,
                logical=logical[order], u=u[order], bcc=bcc[order])


def check(directory, dimension):
    history = athena_read.hst(str(directory / f"{BASENAME}.mhd.hst"))
    user = athena_read.hst(str(directory / f"{BASENAME}.user.hst"))
    for name in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd",
                 "lf_mirror", "lf_firehs", "lf_hwproj"):
        assert np.all(history[name] == 0.0), name
    for name in user:
        if name.startswith("amr_") or name == "bad_state":
            assert np.all(user[name] == 0.0), name
    assert np.max(user["max_ndiv"]) < 1.0e-12
    assert history["lf_nstage"][-1] > 0
    np.testing.assert_allclose(history["tot-E"], history["tot-E"][0],
                               rtol=0.0, atol=5.0e-12)
    paths = sorted((directory / "rst").glob(f"{BASENAME}.*.rst"))
    assert len(paths) >= 2
    states = [read_state(path) for path in paths]
    initial, final = states[0], states[-1]
    assert initial["cycle"] == 0 and final["cycle"] > 0
    assert final["time"] == 0.002
    assert initial["nmb"] == (25 if dimension == 2 else 148)
    for state in states:
        for name in ("energy", "mu"):
            np.testing.assert_allclose(state["totals"][name], initial["totals"][name],
                                       rtol=0.0, atol=5.0e-12, err_msg=name)
        np.testing.assert_array_equal(state["u"][:, [0, 1, 2, 3, 6]],
                                      initial["u"][:, [0, 1, 2, 3, 6]])
        np.testing.assert_array_equal(state["bcc"], initial["bcc"])
    assert np.max(np.abs(final["u"][:, 4] - initial["u"][:, 4])) > 1.0e-8
    return final


def assert_same(reference, candidate):
    assert reference["time"] == candidate["time"]
    assert reference["cycle"] == candidate["cycle"]
    for name in ("logical", "u", "bcc"):
        np.testing.assert_array_equal(reference[name], candidate[name])
