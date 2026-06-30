from __future__ import annotations

import os
from pathlib import Path
import struct
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
LOAD_BALANCE_SOURCE = REPO_ROOT / "src/mesh/load_balance.cpp"


def _load_balance_definition() -> str:
    source = LOAD_BALANCE_SOURCE.read_text(encoding="ascii")
    signature = (
        "void Mesh::LoadBalance(float *clist, int *rlist, int *slist, "
        "int *nlist, int nb)"
    )
    start = source.index(signature)
    brace = source.index("{", start)
    depth = 0
    for index in range(brace, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[start:index + 1]
    raise AssertionError("unterminated Mesh::LoadBalance definition")


@pytest.fixture(scope="session")
def load_balance_harness(tmp_path_factory: pytest.TempPathFactory) -> Path:
    build_dir = tmp_path_factory.mktemp("mesh_load_balance")
    harness = build_dir / "mesh_load_balance_harness.cpp"
    executable = build_dir / "mesh_load_balance_harness"
    harness.write_text(
        r"""
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#define MPI_PARALLEL_ENABLED 1

namespace global_variable {
int my_rank = 0;
int nranks = 1;
} // namespace global_variable

class Mesh {
 public:
  bool adaptive = false;
  void LoadBalance(float *clist, int *rlist, int *slist, int *nlist, int nb);
};

"""
        + _load_balance_definition()
        + r"""

int main() {
  int nb = 0;
  if (!(std::cin >> global_variable::nranks >> nb)) return 2;

  std::vector<float> costs(nb);
  for (float &cost : costs) {
    if (!(std::cin >> cost)) return 2;
  }
  std::vector<int> ranks(nb, -1);
  std::vector<int> starts(global_variable::nranks, -1);
  std::vector<int> counts(global_variable::nranks, -1);

  Mesh mesh;
  mesh.LoadBalance(costs.data(), ranks.data(), starts.data(), counts.data(), nb);

  std::cout << "ranks";
  for (int rank : ranks) std::cout << ' ' << rank;
  std::cout << '\n' << "starts";
  for (int start : starts) std::cout << ' ' << start;
  std::cout << '\n' << "counts";
  for (int count : counts) std::cout << ' ' << count;
  std::cout << '\n';
  return 0;
}
""",
        encoding="ascii",
    )
    subprocess.run(
        [
            os.environ.get("CXX", "c++"),
            "-std=c++17",
            "-Wall",
            "-Wextra",
            "-Werror",
            str(harness),
            "-o",
            str(executable),
        ],
        cwd=REPO_ROOT,
        check=True,
    )
    return executable


def _invoke(
    executable: Path, costs: list[float], nranks: int
) -> subprocess.CompletedProcess[str]:
    payload = f"{nranks} {len(costs)}\n" + " ".join(
        format(cost, ".9g") for cost in costs
    )
    return subprocess.run(
        [str(executable)],
        cwd=REPO_ROOT,
        input=payload,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def _partition(
    executable: Path, costs: list[float], nranks: int
) -> tuple[list[int], list[int], list[int]]:
    result = _invoke(executable, costs, nranks)
    assert result.returncode == 0, result.stdout + result.stderr
    vectors: dict[str, list[int]] = {}
    for line in result.stdout.splitlines():
        fields = line.split()
        if fields and fields[0] in {"ranks", "starts", "counts"}:
            vectors[fields[0]] = [int(value) for value in fields[1:]]
    assert vectors.keys() == {"ranks", "starts", "counts"}
    return vectors["ranks"], vectors["starts"], vectors["counts"]


def _assert_complete_partition(
    partition: tuple[list[int], list[int], list[int]], nb: int, nranks: int
) -> None:
    ranks, starts, counts = partition
    assert len(ranks) == nb
    assert len(starts) == nranks
    assert len(counts) == nranks
    assert all(count > 0 for count in counts)
    assert sum(counts) == nb
    assert sorted(set(ranks)) == list(range(nranks))

    next_start = 0
    for rank, (start, count) in enumerate(zip(starts, counts)):
        assert start == next_start
        assert ranks[start:start + count] == [rank] * count
        next_start += count
    assert next_start == nb


def _float32(value: float) -> float:
    return struct.unpack("=f", struct.pack("=f", value))[0]


def _legacy_rank_ids(costs: list[float], nranks: int) -> list[int]:
    costs = [_float32(cost) for cost in costs]
    total_cost = _float32(0.0)
    for cost in costs:
        total_cost = _float32(total_cost + cost)

    rank = nranks - 1
    target_cost = _float32(total_cost / nranks)
    rank_cost = _float32(0.0)
    rank_ids = [-1] * len(costs)
    for index in range(len(costs) - 1, -1, -1):
        rank_cost = _float32(rank_cost + costs[index])
        rank_ids[index] = rank
        if rank_cost >= target_cost and rank > 0:
            rank -= 1
            total_cost = _float32(total_cost - rank_cost)
            rank_cost = _float32(0.0)
            target_cost = _float32(total_cost / (rank + 1))
    return rank_ids


@pytest.mark.parametrize(
    "costs",
    [
        pytest.param([100.0, 1.0, 1.0, 1.0], id="dominant-beginning"),
        pytest.param([1.0, 100.0, 1.0, 1.0], id="dominant-middle"),
        pytest.param([1.0, 1.0, 1.0, 100.0], id="dominant-end"),
    ],
)
def test_four_blocks_four_ranks_reserve_one_block_each(
    load_balance_harness: Path, costs: list[float]
) -> None:
    partition = _partition(load_balance_harness, costs, nranks=4)
    _assert_complete_partition(partition, nb=4, nranks=4)
    assert partition[2] == [1, 1, 1, 1]


def test_more_blocks_than_ranks_keeps_complete_contiguous_partition(
    load_balance_harness: Path,
) -> None:
    costs = [100.0, 1.0, 7.0, 2.0, 11.0, 3.0, 1.0, 5.0, 2.0]
    partition = _partition(load_balance_harness, costs, nranks=4)
    _assert_complete_partition(partition, nb=len(costs), nranks=4)


def test_q011_1742_block_736_rank_heterogeneous_replay(
    load_balance_harness: Path,
) -> None:
    tail = [float(1 + ((index * 37 + 11) % 23)) for index in range(1, 1742)]
    dominant_cost = float(sum(tail) * 736 * 4)
    costs = [dominant_cost, *tail]

    legacy_rank_ids = _legacy_rank_ids(costs, nranks=736)
    assert set(legacy_rank_ids) == {735}

    partition = _partition(load_balance_harness, costs, nranks=736)
    _assert_complete_partition(partition, nb=1742, nranks=736)


def test_more_ranks_than_blocks_fails_clearly(load_balance_harness: Path) -> None:
    result = _invoke(load_balance_harness, [1.0, 1.0, 1.0], nranks=4)
    assert result.returncode != 0
    assert "requires at least one MeshBlock per rank" in result.stdout
