"""Exercise the production load balancer with uneven costs and invalid rank counts."""

import os
from pathlib import Path
import subprocess

import pytest


@pytest.fixture(scope="module")
def load_balance_harness(tmp_path_factory):
    path = Path(__file__).resolve().parents[1] / "src/mesh/load_balance.cpp"
    source = path.read_text()
    start = source.index("\nvoid Mesh::LoadBalance(") + 1
    definition = source[start:source.index("\n//----------------", start)]
    directory = tmp_path_factory.mktemp("load_balance")
    harness = directory / "test.cpp"
    executable = directory / "test"
    harness.write_text(
        """
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>
#define MPI_PARALLEL_ENABLED 1
namespace global_variable { int my_rank = 0, nranks = 1; }
class Mesh {
 public:
  bool adaptive = false;
  void LoadBalance(float *, int *, int *, int *, int);
};
""" + definition + r"""
int main() {
  int nb;
  std::cin >> global_variable::nranks >> nb;
  std::vector<float> costs(nb);
  for (auto &cost : costs) std::cin >> cost;
  std::vector<int> ranks(nb, -1), starts(global_variable::nranks, -1),
                   counts(global_variable::nranks, -1);
  Mesh mesh;
  mesh.LoadBalance(costs.data(), ranks.data(), starts.data(), counts.data(), nb);
  for (auto &values : {ranks, starts, counts}) {
    for (int value : values) std::cout << value << ' ';
    std::cout << '\n';
  }
}
"""
    )
    subprocess.run([os.environ.get("CXX", "c++"), "-std=c++17", "-Wall", "-Wextra",
                    "-Werror", str(harness), "-o", str(executable)], check=True)
    return executable


@pytest.mark.parametrize("costs,nranks", [
    ([100.0, 1.0, 1.0, 1.0], 4),
    ([1.0, 100.0, 1.0, 1.0], 4),
    ([1.0, 1.0, 1.0, 100.0], 4),
    ([100.0, 1.0, 7.0, 2.0, 11.0, 3.0, 1.0, 5.0, 2.0], 4),
    ([1.0e9] + [float(1 + ((i * 37 + 11) % 23)) for i in range(1, 1742)], 736),
    ([1.0, 1.0, 1.0], 4),
])
def test_load_balance_partition(load_balance_harness, costs, nranks):
    result = subprocess.run(
        [str(load_balance_harness)], input=f"{nranks} {len(costs)}\n" +
        " ".join(map(str, costs)), text=True, capture_output=True,
    )
    if len(costs) < nranks:
        assert result.returncode != 0
        assert "requires at least one MeshBlock per rank" in result.stdout
        return
    assert result.returncode == 0, result.stdout + result.stderr
    ranks, starts, counts = [list(map(int, line.split()))
                             for line in result.stdout.strip().splitlines()[-3:]]
    assert len(ranks) == len(costs)
    assert len(starts) == len(counts) == nranks
    next_start = 0
    for rank, (start, count) in enumerate(zip(starts, counts)):
        assert count > 0 and start == next_start
        assert ranks[start:start + count] == [rank] * count
        next_start += count
    assert next_start == len(costs)
