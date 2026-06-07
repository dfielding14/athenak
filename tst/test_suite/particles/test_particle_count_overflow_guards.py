"""Source-level regression tests for particle MPI/output overflow guards."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[3]


def _source(relative_path):
    return (ROOT / relative_path).read_text(encoding="utf-8")


def test_boundary_exchange_checks_mpi_counts_and_displacements():
    source = _source("src/bvals/bvals_part.cpp")

    assert "CheckedMpiIntProduct" in source
    assert "CheckedMpiIntAdd" in source
    assert "CheckedSizeToMpiInt" in source
    assert "nrdata*(sends_thisrank[n].nprtcls)" not in source
    assert "nidata*(sends_thisrank[n].nprtcls)" not in source
    assert "(pmy_part->nrdata)*(recvs_thisrank[n].nprtcls)" not in source
    assert "(pmy_part->nidata)*(recvs_thisrank[n].nprtcls)" not in source


def test_amr_alltoallv_checks_counts_and_displacements():
    source = _source("src/particles/particles_lagrangian_mc.cpp")

    assert "CheckedParticleIntProduct" in source
    assert "CheckedParticleIntAdd" in source
    for expression in (
        "send_particles[r]*nrdata",
        "recv_particles[r]*nrdata",
        "send_particles[r]*nidata",
        "recv_particles[r]*nidata",
    ):
        assert expression not in source

    alltoallv_calls = re.findall(r"MPI_Alltoallv\([^;]+;", source, re.DOTALL)
    assert len(alltoallv_calls) == 3
    assert all("count.data()" in call and "disp.data()" in call
               for call in alltoallv_calls)


def test_history_gather_checks_record_products_and_displacements():
    source = _source("src/outputs/prtcl_thermo_history.cpp")

    assert "CheckedHistoryIntProduct" in source
    assert "CheckedHistoryIntAdd" in source
    assert "CheckedHistoryIntFromUInt64" in source
    assert "counts[n]*int_per_record" not in source
    assert "counts[n]*real_per_record" not in source
    assert "tag_displ[n-1] + tag_counts[n-1]" not in source


def test_global_particle_totals_use_uint64():
    mesh_header = _source("src/mesh/mesh.hpp")
    outputs_header = _source("src/outputs/outputs.hpp")

    assert "std::uint64_t nprtcl_total" in mesh_header
    assert outputs_header.count("std::uint64_t npout_total") == 2


def test_particle_vtk_rank_offsets_use_uint64():
    source = _source("src/outputs/vtk_prtcl.cpp")

    assert "std::vector<std::uint64_t> rank_offset" in source
    assert "std::vector<int> rank_offset" not in source
    assert "std::size_t local_position_count" in source


def test_ito_fatal_path_aborts_mpi_and_limiter_uses_wide_index():
    source = _source("src/particles/particles_lagrangian_ito.cpp")

    assert "MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE)" in source
    assert "Kokkos::IndexType<std::int64_t>" in source
    assert "std::int64_t nmkji" in source
