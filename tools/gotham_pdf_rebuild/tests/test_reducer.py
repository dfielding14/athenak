"""Executable-level validation for the standalone GOTHAM PDF reducer."""

from __future__ import annotations

import json
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pytest

from reference import (
    ORIGINAL_PRODUCTS,
    SCIENCE_PRODUCT_COUNT,
    SCIENCE_PRODUCTS,
    SCIENCE_TOTAL_BINS,
    Product,
    parse_header,
    reference_histograms,
)
from synthetic_gotham import SEQUENCE_DEFAULT, TIME_DEFAULT, write_dataset


TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parents[2]
DEFAULT_EXE = (
    REPO_ROOT / "build-gotham-pdf-frontier-serial" / "gotham_pdf_rebuild"
)
ANALYSIS_ROOT = Path(
    os.environ.get(
        "GOTHAM_ANALYSIS_ROOT", "/lustre/orion/ast207/proj-shared/gotham/analysis"
    )
)


def reducer_executable() -> Path:
    path = Path(os.environ.get("GOTHAM_PDF_REBUILD_EXE", str(DEFAULT_EXE))).resolve()
    if not path.is_file():
        pytest.skip(f"reducer executable not found: {path}")
    return path


def mpi_launcher() -> str:
    launcher = os.environ.get("GOTHAM_MPIEXEC") or shutil.which("mpiexec") or shutil.which("mpirun")
    if launcher is None:
        pytest.skip("mpiexec/mpirun is unavailable")
    return launcher


def run_reducer(
    input_dir: Path,
    output_dir: Path,
    products: str,
    *,
    ranks: int = 1,
    max_blocks_per_shard: int = 0,
    only_products: Iterable[str] = (),
) -> subprocess.CompletedProcess:
    command = []
    if ranks > 1:
        command.extend((mpi_launcher(), "-n", str(ranks)))
    command.extend(
        (
            str(reducer_executable()),
            "--input-dir",
            str(input_dir),
            "--sequence",
            SEQUENCE_DEFAULT,
            "--output-dir",
            str(output_dir),
            "--products",
            products,
            "--chunk-blocks",
            "1",
        )
    )
    for product in only_products:
        command.extend(("--only-product", product))
    if max_blocks_per_shard:
        command.extend(("--max-blocks-per-shard", str(max_blocks_per_shard)))
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    return subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=environment,
        timeout=600,
    )


def payload_path(output_dir: Path, product: Product) -> Path:
    return output_dir / product.identifier / f"gotham.{SEQUENCE_DEFAULT}.pdf"


def assert_header_contract(output_dir: Path, product: Product) -> None:
    header = parse_header(output_dir / product.identifier / "gotham.header.pdf")
    assert header["format"] == "dense"
    assert header["distribution"] == "global_rebuild"
    assert int(header["ndim"]) == len(product.axes)
    assert int(header["total_bins"]) == product.total_bins
    assert header["weight"] == (
        product.weight if product.weight in {"mass", "volume"} else "variable"
    )
    if product.weight not in {"mass", "volume"}:
        assert header["weight_variable"] == product.weight

    for dimension, (axis, stride) in enumerate(zip(product.axes, product.strides), start=1):
        assert header[f"variable_{dimension}"] == axis.variable
        assert int(header[f"nbin{dimension}"]) == axis.nbin
        assert float(header[f"bin{dimension}_min"]) == pytest.approx(axis.minimum)
        assert float(header[f"bin{dimension}_max"]) == pytest.approx(axis.maximum)
        assert header[f"scale{dimension}"] == axis.scale
        assert int(header[f"stride{dimension}"]) == stride
        edges = np.fromstring(header[f"bin_edges_{dimension}"], sep=" ")
        assert edges.size == axis.nbin + 1
        assert edges[0] == pytest.approx(axis.minimum)
        assert edges[-1] == pytest.approx(axis.maximum)
        assert np.all(np.diff(edges) > 0.0)


def assert_sparse_reference(
    output_dir: Path,
    products: Iterable[Product],
    reference,
    *,
    rtol: float = 2.0e-12,
    atol: float = 1.0e-12,
) -> None:
    for product in products:
        # AthenaK's cooling interpolation uses fused multiply-add operations.
        # The independent Python reference intentionally does not depend on a
        # platform-specific fma binding, so luminosity-weighted cooling sums
        # need a slightly wider relative tolerance than the other products.
        product_rtol = 2.0e-8 if product.weight == "edot_cool" else rtol
        path = payload_path(output_dir, product)
        raw = np.memmap(path, dtype=np.float64, mode="r")
        assert raw.size == product.total_bins + 1
        assert raw[0] == pytest.approx(TIME_DEFAULT)
        actual = raw[1:]
        expected_sparse = reference[product.identifier]
        expected_indices = np.fromiter(expected_sparse, dtype=np.int64)
        expected_values = np.fromiter(expected_sparse.values(), dtype=np.float64)
        expected = np.zeros(product.total_bins, dtype=np.float64)
        expected[expected_indices] = expected_values
        assert np.allclose(actual, expected, rtol=product_rtol, atol=atol)
        assert np.isclose(
            actual.sum(dtype=np.float64),
            expected_values.sum(dtype=np.float64),
            rtol=product_rtol,
            atol=atol,
        )
        del expected
        del raw


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    root = tmp_path_factory.mktemp("gotham-synthetic")
    blocks = write_dataset(root / "input", blocks=4, shards=2)
    return root, blocks


@pytest.fixture(scope="module")
def science_serial(synthetic):
    root, blocks = synthetic
    serial = root / "science-serial"
    selected = [product.identifier for product in SCIENCE_PRODUCTS]
    run_reducer(
        root / "input", serial, "science", ranks=1, only_products=selected
    )
    return serial, blocks


@pytest.fixture(scope="module")
def science_mpi2(synthetic):
    root, _ = synthetic
    mpi2 = root / "science-mpi2"
    selected = [product.identifier for product in SCIENCE_PRODUCTS]
    run_reducer(
        root / "input", mpi2, "science", ranks=2, only_products=selected
    )
    return mpi2


def test_science_products_match_independent_reference(science_serial) -> None:
    serial, blocks = science_serial
    reference = reference_histograms(blocks, SCIENCE_PRODUCTS)
    for product in SCIENCE_PRODUCTS:
        assert_header_contract(serial, product)
    assert_sparse_reference(serial, SCIENCE_PRODUCTS, reference)

    manifest = json.loads((serial / "rebuild_manifest.json").read_text(encoding="utf-8"))
    assert manifest["shards_available"] == 2
    assert manifest["shards_processed"] == 2
    assert manifest["meshblocks_processed"] == len(blocks)
    assert manifest["cells_processed"] == sum(block.cells for block in blocks)


def test_science_products_are_mpi_invariant(science_serial, science_mpi2) -> None:
    serial, _ = science_serial
    mpi2 = science_mpi2
    for product in SCIENCE_PRODUCTS:
        assert (serial / product.identifier / "gotham.header.pdf").read_bytes() == (
            mpi2 / product.identifier / "gotham.header.pdf"
        ).read_bytes()
        serial_values = np.memmap(payload_path(serial, product), dtype=np.float64, mode="r")
        mpi_values = np.memmap(payload_path(mpi2, product), dtype=np.float64, mode="r")
        assert np.allclose(serial_values, mpi_values, rtol=2.0e-12, atol=1.0e-12)
        del serial_values
        del mpi_values


def test_underflow_overflow_and_signed_weights_are_exercised(science_serial) -> None:
    serial, _ = science_serial
    reader_path = payload_path(serial, SCIENCE_PRODUCTS[0])
    values = np.memmap(reader_path, dtype=np.float64, mode="r", offset=8).reshape(
        SCIENCE_PRODUCTS[0].shape
    )
    assert np.count_nonzero(values[0, :, :]) == 0
    assert np.count_nonzero(values[-1, :, :]) > 0
    assert np.count_nonzero(values[:, 0, :]) == 0
    assert np.count_nonzero(values[:, -1, :]) == 0
    assert np.count_nonzero(values[:, :, 0]) > 0
    assert np.count_nonzero(values[:, :, -1]) > 0
    del values

    inflow_energy = np.memmap(
        payload_path(
            serial,
            next(
                p
                for p in SCIENCE_PRODUCTS
                if p.identifier.endswith("edot_in_abs")
            ),
        ),
        dtype=np.float64,
        mode="r",
        offset=8,
    )
    vertical_ram = np.memmap(
        payload_path(
            serial,
            next(
                p
                for p in SCIENCE_PRODUCTS
                if p.identifier.endswith("vertical_geometry_ram_out")
            ),
        ),
        dtype=np.float64,
        mode="r",
        offset=8,
    )
    assert np.all(inflow_energy >= 0.0)
    assert np.all(vertical_ram >= 0.0)
    assert inflow_energy.sum() > 0.0
    assert vertical_ram.sum() > 0.0
    del inflow_energy
    del vertical_ram


def test_science_catalog_count_and_memory_contract() -> None:
    completed = subprocess.run(
        [str(reducer_executable()), "--products", "science", "--list-products"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=60,
    )
    product_lines = [
        line
        for line in completed.stdout.splitlines()
        if line.startswith("science_")
    ]
    assert len(product_lines) == SCIENCE_PRODUCT_COUNT
    assert f"total_bins={SCIENCE_TOTAL_BINS}" in completed.stdout
    assert not any("output" in line for line in product_lines)


def test_runtime_64_cubed_meshblock_is_processed(tmp_path: Path) -> None:
    input_dir = tmp_path / "input64"
    blocks = write_dataset(input_dir, blocks=1, shards=1, nx1=64)
    output = tmp_path / "science64"
    representative = SCIENCE_PRODUCTS[0]
    run_reducer(
        input_dir,
        output,
        "science",
        ranks=1,
        only_products=[representative.identifier],
    )

    manifest = json.loads((output / "rebuild_manifest.json").read_text(encoding="utf-8"))
    assert manifest["meshblocks_processed"] == 1
    assert manifest["cells_processed"] == 64**3
    reference = reference_histograms(blocks, [representative])
    assert_sparse_reference(output, [representative], reference)


def test_geometry_must_match_logical_key(tmp_path: Path) -> None:
    input_dir = tmp_path / "invalid-geometry"
    write_dataset(input_dir, blocks=1, shards=1)
    shard = input_dir / "node_00000000" / f"gotham.hydro_w.{SEQUENCE_DEFAULT}.bin"
    with shard.open("r+b") as stream:
        while True:
            line = stream.readline()
            if line.startswith(b"  header offset="):
                embedded_bytes = int(line.split(b"=", 1)[1])
                break
        stream.seek(embedded_bytes, os.SEEK_CUR)
        geometry_offset = stream.tell() + 10 * 4
        stream.seek(geometry_offset)
        x1min = struct.unpack("=d", stream.read(8))[0]
        stream.seek(geometry_offset)
        stream.write(struct.pack("=d", x1min + 1.0))

    with pytest.raises(subprocess.CalledProcessError):
        run_reducer(input_dir, tmp_path / "should-fail", "original")


@pytest.mark.large
def test_every_original_stream_and_analysis_reader(synthetic) -> None:
    root, blocks = synthetic
    output = root / "original"
    run_reducer(root / "input", output, "original", ranks=1)
    reference = reference_histograms(blocks, ORIGINAL_PRODUCTS)
    for product in ORIGINAL_PRODUCTS:
        assert_header_contract(output, product)
    assert_sparse_reference(output, ORIGINAL_PRODUCTS, reference)

    if not (ANALYSIS_ROOT / "gotham_analysis" / "readers" / "pdf.py").is_file():
        pytest.skip(f"production analysis reader unavailable beneath {ANALYSIS_ROOT}")
    sys.path.insert(0, str(ANALYSIS_ROOT))
    try:
        from gotham_analysis.readers.pdf import read_pdf

        representative = next(p for p in ORIGINAL_PRODUCTS if p.identifier == "output29")
        loaded = read_pdf(payload_path(output, representative))
        assert loaded["time"] == pytest.approx(TIME_DEFAULT)
        assert loaded["header"]["shape"] == representative.shape
        assert loaded["header"]["variables"] == tuple(axis.variable for axis in representative.axes)
        assert loaded["pdf"].shape == representative.shape
    finally:
        sys.path.remove(str(ANALYSIS_ROOT))
