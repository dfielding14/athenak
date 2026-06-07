"""Helpers for constructing historical Ito particle restart layouts."""

import struct


PARTICLE_HEADER = struct.Struct("@16s6iQ")


def make_legacy_restart(source, destination, version):
    """Convert a current Ito restart into the historical v2 or v3 layout."""
    payload = bytearray(source.read_bytes())
    marker = b"<par_end>"
    parameter_end = payload.index(marker) + len(marker) + 1
    parameter_lines = payload[:parameter_end].decode("utf-8").splitlines(keepends=True)
    parameter_lines = [
        line for line in parameter_lines
        if not line.lstrip().startswith("ito_covariance_model")
    ]
    payload = (
        bytearray("".join(parameter_lines).encode("utf-8"))
        + payload[parameter_end:]
    )

    particle_offset = payload.rfind(b"ATHKPRTCLMC")
    fields = PARTICLE_HEADER.unpack_from(payload, particle_offset)
    assert fields[1] == 4
    assert fields[2] == 2
    _, _, _, nrdata, nidata, nlocal, nschedules, _ = fields
    int_size = struct.calcsize("@i")
    tag_size = struct.calcsize("@Q")
    fixed_bytes = (
        PARTICLE_HEADER.size
        + int_size
        + 3 * nschedules * int_size
        + nidata * nlocal * int_size
        + nlocal * tag_size
    )
    real_count = nschedules + nrdata * nlocal
    real_size, remainder = divmod(
        len(payload) - particle_offset - fixed_bytes, real_count
    )
    assert remainder == 0
    assert real_size in (4, 8)
    real_base = (
        particle_offset
        + PARTICLE_HEADER.size
        + int_size
        + nschedules * real_size
        + 3 * nschedules * int_size
    )
    int_base = real_base + nrdata * nlocal * real_size
    tag_base = int_base + nidata * nlocal * int_size
    assert tag_base + nlocal * tag_size == len(payload)
    if version == 2:
        for particle in range(nlocal):
            tag = struct.unpack_from("@Q", payload, tag_base + particle * tag_size)[0]
            assert tag <= 2**31 - 1
            struct.pack_into(
                "@i", payload, int_base + (nlocal + particle) * int_size, tag
            )
        del payload[tag_base:tag_base + nlocal * tag_size]
    del payload[
        particle_offset + PARTICLE_HEADER.size:
        particle_offset + PARTICLE_HEADER.size + int_size
    ]
    struct.pack_into("@i", payload, particle_offset + 16, version)
    destination.write_bytes(payload)
