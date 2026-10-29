"""Exact cell/face redistribution and its composition with uniform refinement."""

import re
import shutil

import numpy as np
import pytest

from test_suite.cgl.test_refine_on_restart import (
    MODES, global_fields, input_text, mpi_available, patch_checkpoint, run,
)


def block_flags(shape):
    return tuple(f"meshblock/nx{axis}={cells}" for axis, cells in enumerate(shape, 1))


def assert_same_fields(first, second):
    assert first["shape"] == second["shape"]
    first_u, first_b = global_fields(first)
    second_u, second_b = global_fields(second)
    np.testing.assert_array_equal(first_u, second_u)
    for a, b in zip(first_b, second_b):
        np.testing.assert_array_equal(a, b)
    assert (first["time"], first["cycle"]) == (second["time"], second["cycle"])


@pytest.mark.parametrize("mode", MODES)
def test_reblock_exact_split_merge_and_round_trip(tmp_path, mode):
    mpi = mpi_available()
    initial = run(tmp_path / "initial", mode,
                  source=input_text(mode, blocks=2, nghost=2),
                  flags=("mhd/fofc=false",))
    clean = patch_checkpoint(initial, mode, tmp_path / "clean.rst", nonuniform=True)
    dirty = patch_checkpoint(initial, mode, tmp_path / "dirty.rst",
                             nonuniform=True, poison_ghosts=True)
    reference = run(tmp_path / "reference", mode, restart=clean)
    split = run(tmp_path / "split", mode, restart=dirty, reblock=True,
                flags=block_flags((4, 4, 4)), ranks=4 if mpi else 1)
    assert split["block_shape"] == (4, 4, 4) and len(split["u"]) == 16
    assert_same_fields(split, reference)
    assert re.search(r"^reblock_on_restart\s*=\s*(?:false|0)\s",
                     split["text"], re.MULTILINE)
    held = run(tmp_path / "held", mode, restart=split["path"], ranks=2 if mpi else 1)
    assert held["block_shape"] == split["block_shape"]
    assert_same_fields(held, split)
    merged = run(tmp_path / "merged", mode, restart=held["path"], reblock=True,
                 flags=block_flags((16, 8, 8)))
    assert merged["block_shape"] == (16, 8, 8) and len(merged["u"]) == 1
    assert_same_fields(merged, reference)
    restored = run(tmp_path / "restored", mode, restart=merged["path"], reblock=True,
                   flags=block_flags((8, 8, 8)), ranks=2 if mpi else 1)
    assert_same_fields(restored, reference)


def test_reblock_exact_non_nested_blocks_and_non_cubic_cells(tmp_path):
    mpi = mpi_available()
    source = input_text("ideal", shape=(30, 10, 10), block_shape=(6, 5, 5))
    source += "\n<mesh>\nx1max = 1\n"
    initial = run(tmp_path / "initial", "ideal", source=source,
                  flags=("mhd/fofc=false",))
    clean = patch_checkpoint(initial, "ideal", tmp_path / "clean.rst", nonuniform=True)
    dirty = patch_checkpoint(initial, "ideal", tmp_path / "dirty.rst",
                             nonuniform=True, poison_ghosts=True)
    reference = run(tmp_path / "reference", "ideal", restart=clean)
    result = run(tmp_path / "reblocked", "ideal", restart=dirty, reblock=True,
                 flags=block_flags((5, 10, 5)), ranks=3 if mpi else 1)
    assert result["block_shape"] == (5, 10, 5) and len(result["u"]) == 12
    assert_same_fields(result, reference)
    returned = run(tmp_path / "returned", "ideal", restart=result["path"],
                   reblock=True, flags=block_flags((6, 5, 5)), ranks=2 if mpi else 1)
    assert_same_fields(returned, reference)


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("block_shape", [(32, 16, 16), (4, 4, 4)])
def test_refine_and_reblock_matches_fixed_block_refinement(tmp_path, mode, block_shape):
    mpi = mpi_available()
    initial = run(tmp_path / "initial", mode, source=input_text(mode, blocks=2))
    source = patch_checkpoint(initial, mode, tmp_path / "source.rst", nonuniform=True)
    fixed = run(tmp_path / "fixed", mode, restart=source, refine=True,
                ranks=2 if mpi else 1)
    combined = run(tmp_path / "combined", mode, restart=source,
                   refine=True, reblock=True, flags=block_flags(block_shape),
                   ranks=2 if mpi and block_shape == (4, 4, 4) else 1)
    assert combined["block_shape"] == block_shape
    assert len(combined["u"]) == np.prod(combined["shape"])//np.prod(block_shape)
    expected_u, expected_b = global_fields(fixed)
    actual_u, actual_b = global_fields(combined)
    np.testing.assert_allclose(actual_u, expected_u, rtol=5e-13, atol=5e-13)
    for actual, expected in zip(actual_b, expected_b):
        np.testing.assert_array_equal(actual, expected)
    assert (combined["time"], combined["cycle"]) == (fixed["time"], fixed["cycle"])
    assert re.search(r"^refine_on_restart\s*=\s*(?:false|0)\s",
                     combined["text"], re.MULTILINE)


@pytest.mark.parametrize("flags,message", [
    (block_flags((6, 8, 8)), "evenly divisible"),
    (block_flags((2, 8, 8)), "MeshBlock must be >= 4"),
    (("mesh/nghost=2",), "ghost"),
    (("mesh/nx1=16",), "mesh dimensions unchanged"),
])
def test_reblock_rejects_incompatible_geometry(tmp_path, flags, message):
    initial = run(tmp_path / "initial", "ideal", source=input_text("ideal"))
    output = run(tmp_path / "invalid", "ideal", restart=initial["path"],
                 reblock=True, flags=flags, error=True)
    assert message in output


def test_reblock_requires_shared_restart(tmp_path):
    output = run(tmp_path / "fresh", "ideal", source=input_text("ideal"),
                 reblock=True, error=True)
    assert "requires a shared checkpoint" in output
    initial = run(tmp_path / "initial", "ideal", source=input_text("ideal"))
    directory = tmp_path / "rank_00000000"
    directory.mkdir()
    checkpoint = directory / initial["path"].name
    shutil.copyfile(initial["path"], checkpoint)
    output = run(tmp_path / "rank_local", "ideal", restart=checkpoint,
                 reblock=True, ranks=2, error=True)
    assert "requires a shared checkpoint" in output


def test_reblock_checks_final_rank_count(tmp_path):
    initial = run(tmp_path / "initial", "ideal", source=input_text("ideal", blocks=2))
    output = run(tmp_path / "invalid", "ideal", restart=initial["path"], reblock=True,
                 flags=block_flags((16, 8, 8)), ranks=2, error=True)
    assert "Fewer MeshBlocks (nmb_total=1) than MPI ranks (nranks=2)" in output


def test_reblock_rejects_amr_checkpoint(tmp_path):
    source = input_text("cgl") + """
<mesh_refinement>
refinement = adaptive
num_levels = 2
prolong_primitives = true
max_nmb_per_rank = 8
<amr_criterion1>
method = time
time_start = 1
"""
    initial = run(tmp_path / "initial", "cgl", source=source)
    output = run(tmp_path / "invalid", "cgl", restart=initial["path"], reblock=True,
                 flags=block_flags((4, 4, 4)), error=True)
    assert "without SMR/AMR" in output
