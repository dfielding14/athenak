"""Exercise the production scalar closure independently of the MHD integrator."""

from pathlib import Path
import subprocess

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[3]


@pytest.fixture(scope="module")
def closure(tmp_path_factory):
    directory = tmp_path_factory.mktemp("current_limited_closure")
    # The unused edge kernels only need these types declared; no mock closure.
    (directory / "athena.hpp").write_text("""
#pragma once
#include <array>
#include <initializer_list>
#include <stdexcept>
#define KOKKOS_INLINE_FUNCTION inline
using Real = double;
constexpr int IDN = 0;
template <typename T> struct Array {
  mutable std::array<T,4096> data{};
  int extent_int(int) const { return 0; }
  template <typename... I> T &operator()(I... index) const {
    int offset=0;
    for (int i : {index...}) {
      if (i<0 || i>=8) throw std::out_of_range("stencil exceeds two ghost cells");
      offset=8*offset+i;
    }
    return data.at(offset);
  }
};
template <typename T> using DvceArray5D = Array<T>;
template <typename T> using DvceArray4D = Array<T>;
template <typename T> struct DvceFaceFld4D { Array<T> x1f, x2f, x3f; };
""")
    (directory / "mesh").mkdir()
    (directory / "mesh/mesh.hpp").write_text(
        "#pragma once\nstruct RegionSize { double dx1, dx2, dx3; };\n")
    source = directory / "closure.cpp"
    source.write_text("""
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "diffusion/current_limited_resistivity.hpp"
int main(int argc, char **argv) {
  if (std::string(argv[1]) == "stencil") {
    int dim=std::stoi(argv[2]);
    RegionSize size{0.2,0.3,0.4};
    DvceFaceFld4D<Real> b;
    DvceArray5D<Real> w;
    for (int k=0;k<8;++k) for (int j=0;j<8;++j) for (int i=0;i<8;++i) {
      double x=(i+0.5)*size.dx1,y=(j+0.5)*size.dx2,z=(k+0.5)*size.dx3;
      b.x1f(0,k,j,i)=2*z-3*y;
      b.x2f(0,k,j,i)=5*x-7*z;
      b.x3f(0,k,j,i)=11*y-13*x;
      w(0,IDN,k,j,i)=-1;
    }
    current_limited::Parameters p{0.1,1,1,1,1-sqrt(0.1),sqrt(0.1)};
    std::cout << std::setprecision(17);
    const int js=dim>1?2:0,je=dim>1?5:0,ks=dim>2?2:0,ke=dim>2?5:0;
    for (int c=0;c<3;++c) for (int k=ks;k<=ke+(c!=2);++k)
    for (int j=js;j<=je+(c!=1);++j) for (int i=2;i<=5+(c!=0);++i) {
      auto s=current_limited::EdgeState(b,w,size,p,4,dim>1,dim>2,c,0,k,j,i);
      std::cout << s.j1 << " " << s.j2 << " " << s.j3 << " " << s.q << "\\n";
    }
    for (int k=ks;k<=ke;++k) for (int j=js;j<=je;++j) for (int i=2;i<=5;++i) {
      auto s=current_limited::CellState(b,w,size,p,4,dim>1,dim>2,0,k,j,i);
      std::cout << s.j1 << " " << s.j2 << " " << s.j3 << " " << s.q << "\\n";
    }
    return 0;
  }
  double eta0=std::stod(argv[1]), eta_max=std::stod(argv[2]);
  current_limited::Parameters p{eta0,eta_max,1,1,
      1-sqrt(eta0/eta_max),sqrt(eta0)*sqrt(eta_max)};
  std::cout << std::setprecision(17);
  std::string line;
  while (std::getline(std::cin,line)) {
    std::cout << current_limited::Eta(std::strtod(line.c_str(),nullptr),p) << "\\n";
  }
}
""")
    binary = directory / "closure"
    subprocess.run(["c++", "-std=c++17", "-O2", f"-I{directory}", f"-I{ROOT / 'src'}",
                    str(source), "-o", str(binary)], check=True)
    return binary


def evaluate(binary, q, eta0, eta_max=1):
    result = subprocess.run([str(binary), str(eta0), str(eta_max)],
                            input="\n".join(map(str, q)), text=True,
                            capture_output=True, check=True)
    return np.fromstring(result.stdout, sep="\n")


@pytest.mark.parametrize("ratio", [1e-8, 1e-4, 1e-2, 1])
def test_dense_closure_and_differential_bound(closure, ratio):
    qs = 1 - np.sqrt(ratio)
    q = np.unique(np.r_[0, np.linspace(0, 1, 4001), np.geomspace(1e-12, 1e8, 4001),
                        qs, np.nextafter(qs, 0), np.nextafter(qs, np.inf)])
    eta = evaluate(closure, q, ratio)
    assert eta.shape == q.shape
    assert np.isfinite(eta).all()
    assert eta[0] == ratio
    assert np.min(np.diff(eta)) >= -2e-15
    assert eta.min() >= ratio * (1 - 1e-14)
    assert eta.max() <= 1
    # Secants of E(q)=eta(q)q cannot exceed its maximum derivative. Exclude
    # nextafter-sized gaps, where subtraction cannot resolve the derivative.
    gap = np.diff(q)
    resolved = gap > 1e-3 * np.maximum(q[1:], 1e-12)
    slope = np.diff(q * eta)[resolved] / gap[resolved]
    assert slope.min() >= 0
    assert slope.max() <= 1 + 1e-12
    at_transition = evaluate(closure, [qs], ratio)[0]
    assert at_transition == pytest.approx(np.sqrt(ratio), rel=2e-12)
    if ratio < 1:
        h = (1 - qs) * 1e-3
        y = evaluate(closure, [qs - 2*h, qs - h, qs + h, qs + 2*h], ratio)
        np.testing.assert_allclose([(y[1]-y[0])/h, (y[3]-y[2])/h], 1,
                                   rtol=4e-3, atol=1e-8)


def test_extreme_finite_inputs_remain_finite(closure):
    q = [0, 1e-100, 0.5, 1, 1e100, 1e300]
    for eta0, eta_max in [(1e-200, 1e100), (1e100, 1e200), (1e-200, 1e-200)]:
        eta = evaluate(closure, q, eta0, eta_max)
        assert np.isfinite(eta).all()
        assert (eta >= eta0).all()
        assert (eta <= eta_max).all()


@pytest.mark.parametrize("dim, expected", [
    (1, [0, 13, 5]), (2, [11, 13, 8]), (3, [18, 15, 8]),
])
def test_edge_current_interpolation_and_density_floor(closure, dim, expected):
    result = subprocess.run([str(closure), "stencil", str(dim)], text=True,
                            capture_output=True, check=True)
    rows = np.fromstring(result.stdout, sep=" ").reshape(-1, 4)
    assert len(rows) == {1: 40, 2: 121, 3: 364}[dim]
    np.testing.assert_allclose(rows[:, :3], np.tile(expected, (len(rows), 1)),
                               atol=5e-14, rtol=0)
    np.testing.assert_allclose(rows[:, 3], np.linalg.norm(expected) / 2, rtol=2e-15)
