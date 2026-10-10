"""Execute both production budget loops with float/double process timesteps."""

import os
from pathlib import Path
import shlex
import subprocess

import pytest


@pytest.mark.parametrize("precision", ["float", "double"])
@pytest.mark.parametrize("inactive", [False, True])
def test_parabolic_budget_scaling(tmp_path, precision, inactive):
    root = Path(__file__).resolve().parents[3]
    source = (root / "src/mesh/mesh.cpp").read_text()
    initial = source[source.index("  // Parabolic-process timestep budgets"):
                     source.index("  // z4c timestep")]
    refresh = source[source.index("void Mesh::RefreshSTSParabolicTimeStep() {"):]
    refresh = refresh[refresh.index("  dt_parabolic_sts ="):refresh.index("#if")]
    (tmp_path / "athena.hpp").write_text(f"using Real = {precision};\n")
    unit = tmp_path / "budget.cpp"
    unit.write_text("""
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include "diffusion/parabolic_process.hpp"
using namespace parabolic;
struct Pack { std::vector<ParabolicProcessDescriptor> parabolic_processes; };
int main() {
  const Real sentinel = std::numeric_limits<float>::max();
  for (Real safety : {Real(.3), Real(.9), Real(1)}) {
    const Real cfl_no = .2, sts_safety = safety;
    for (Real raw : {sentinel, std::numeric_limits<Real>::max(),
                     std::numeric_limits<Real>::infinity()}) {
      Real active = .25, explicit_dt = .5;
      Pack pack;
      pack.parabolic_processes.push_back({"inactive", ParabolicProcessOwner::mhd,
          ParabolicIntegratorMode::sts, ParabolicUpdateShape::cell_centered, &raw});
      if (!INACTIVE) {
        pack.parabolic_processes.push_back({"active", ParabolicProcessOwner::mhd,
            ParabolicIntegratorMode::sts, ParabolicUpdateShape::cell_centered, &active});
        pack.parabolic_processes.push_back({"explicit", ParabolicProcessOwner::hydro,
            ParabolicIntegratorMode::explicit_mode,
            ParabolicUpdateShape::cell_centered, &explicit_dt});
      }
      Pack *pmb_pack = &pack;
      Real dt_legacy = sentinel, dt_cycle_candidate = sentinel;
      Real dt_parabolic_sts = sentinel;
""".replace("INACTIVE", "true" if inactive else "false") + initial + """
      const Real expected_sts = INACTIVE ? sentinel : safety*active;
      const Real expected_explicit = INACTIVE ? sentinel : cfl_no*explicit_dt;
      const Real expected_legacy = INACTIVE ? sentinel : cfl_no*active;
      if (dt_parabolic_sts != expected_sts || dt_cycle_candidate != expected_explicit ||
          dt_legacy != expected_legacy) {
        std::cerr << "initial budget lost inactive sentinel or independent safety\\n";
        return 1;
      }
""".replace("INACTIVE", "true" if inactive else "false") + refresh + """
      if (dt_parabolic_sts != expected_sts) {
        std::cerr << "refreshed budget lost inactive sentinel or independent safety\\n";
        return 2;
      }
    }
  }
}
""")
    executable = tmp_path / "budget"
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    subprocess.run(compiler + ["-std=c++17", "-O2", f"-I{tmp_path}",
                              f"-I{root / 'src'}", str(unit), "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    result = subprocess.run([str(executable)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
