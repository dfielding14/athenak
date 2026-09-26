"""Static geometry checks of the production jump estimator on synthetic bcc arrays."""

import math
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[3]


@pytest.fixture(scope="module")
def estimator(tmp_path_factory):
    directory = tmp_path_factory.mktemp("current_limited_estimator")
    # Sparse, checked CPU arrays hold only the populated stencil lines. An access
    # in an inactive direction or outside those lines fails instead of reading zero.
    (directory / "athena.hpp").write_text(r"""
#pragma once
#include <array>
#include <map>
#include <memory>
#define KOKKOS_INLINE_FUNCTION inline
using Real = double;
constexpr int IDN = 0;
template <typename T> struct Array {
  using Key = std::array<int,5>;
  std::shared_ptr<std::map<Key,T>> values = std::make_shared<std::map<Key,T>>();
  T &operator()(int m,int n,int k,int j,int i) const {
    return values->at({m,n,k,j,i});
  }
  T &operator()(int m,int k,int j,int i) const {
    return values->at({m,0,k,j,i});
  }
  void set(int n,int k,int j,int i,T value) { (*values)[{0,n,k,j,i}]=value; }
  void set4(int m,int k,int j,int i,T value) { (*values)[{m,0,k,j,i}]=value; }
  int extent_int(int) const { return values->empty() ? 0 : 65; }
  size_t extent(int n) const { return extent_int(n); }
  size_t size() const { return values->size(); }
};
template <typename T> using DvceArray4D = Array<T>;
template <typename T> using DvceArray5D = Array<T>;
template <typename T> struct DvceFaceFld4D { Array<T> x1f, x2f, x3f; };
""")
    (directory / "mesh").mkdir()
    (directory / "mesh/mesh.hpp").write_text(
        "#pragma once\nstruct RegionSize { double dx1, dx2, dx3; };\n")
    source = directory / "estimator.cpp"
    source.write_text(r"""
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include "diffusion/current_limited_resistivity.hpp"
int main(int argc,char **argv) {
  if (std::string(argv[1])=="cache") {
    const int dim=std::stoi(argv[2]), ng=std::stoi(argv[3]);
    const int ny=dim>1?6:1, nz=dim>2?6:1;
    current_limited::Parameters p{};
    p.b_rec_offset=ng-1;
    for (int m=0;m<2;++m) for (int k=0;k<nz;++k)
    for (int j=0;j<ny;++j) for (int i=0;i<6;++i) {
      p.b_rec_cells.set4(m,k,j,i,10+100*m+i+3*j+7*k);
    }
    double error=0;
    for (int m=0;m<2;++m) for (int k=0;k<nz;++k)
    for (int j=0;j<ny;++j) for (int i=0;i<6;++i) {
      const double actual=current_limited::CellBRec(p,dim>1,dim>2,m,
          k+(dim>2?ng-1:0),j+(dim>1?ng-1:0),i+ng-1);
      error=std::max(error,std::abs(actual-(10+100*m+i+3*j+7*k)));
    }
    current_limited::Parameters constant{};
    constant.b_rec=7;
    std::cout << error << " " << current_limited::CellBRec(
        constant,dim>1,dim>2,1,100,100,100) << "\n";
    return 0;
  }
  std::string profile(argv[1]);
  const int dim=std::stoi(argv[2]);
  const double angle=std::stod(argv[3])*M_PI/180;
  const double shear=std::stod(argv[4])*M_PI/180;
  const double guide=std::stod(argv[5]);
  const double radius=std::stod(argv[6]), floor=std::stod(argv[7]);
  const double amp=std::stod(argv[8]);
  const double dx[3]={std::stod(argv[9]),std::stod(argv[10]),std::stod(argv[11])};
  RegionSize size{dx[0],dx[1],dx[2]};
  double normal[3]={cos(angle),sin(angle),0};
  if (dim==3) for (int d=0;d<3;++d) normal[d]=1/sqrt(3.0);
  const double nxy=hypot(normal[0],normal[1]);
  const double tangent[3]={-normal[1]/nxy,normal[0]/nxy,0};
  const double guide_dir[3]={-normal[2]*tangent[1],normal[2]*tangent[0],
                            normal[0]*tangent[1]-normal[1]*tangent[0]};
  const int center[3]={32,dim>1?32:0,dim>2?32:0};
  DvceArray5D<Real> bcc;
  for (int axis=0;axis<dim;++axis) for (int offset=-32;offset<=32;++offset) {
    int index[3]={center[0],center[1],center[2]};
    index[axis]+=offset;
    const double distance=normal[axis]*offset*dx[axis];
    const double rotation=0.5*shear*tanh(distance);
    for (int c=0;c<3;++c) {
      double value=guide*guide_dir[c];
      if (profile=="harris") value+=amp*tanh(distance)*tangent[c];
      else if (profile=="rotation") {
        value+=amp*(sin(rotation)*tangent[c]+cos(rotation)*guide_dir[c]);
      } else if (profile=="noise") {
        value+=amp*sin(0.5*M_PI*(index[0]+index[1]+index[2]))*tangent[c];
      } else return 2;
      bcc.set(c,index[2],index[1],index[0],value);
    }
  }
  std::cout << std::setprecision(17) << current_limited::JumpEstimate(
      bcc,size,radius,floor,dim>1,dim>2,0,center[2],center[1],center[0]) << "\n";
}
""")
    binary = directory / "estimator"
    subprocess.run(["c++", "-std=c++17", "-O2", f"-I{directory}", f"-I{ROOT / 'src'}",
                    str(source), "-o", str(binary)], check=True)
    return binary


def measure(estimator, profile="harris", dim=1, angle=0, shear=180, guide=0,
            radius=3, floor=0.05, amp=1, spacing=(0.1, 0.1, 0.1)):
    result = subprocess.run([str(estimator), profile, str(dim), str(angle), str(shear),
                             str(guide), str(radius), str(floor), str(amp),
                             *map(str, spacing)], capture_output=True, text=True,
                            check=True)
    return float(result.stdout)


@pytest.mark.parametrize("guide", [0, 1])
def test_harris_center_and_uniform_guide(estimator, guide):
    measured = measure(estimator, guide=guide)
    assert measured == pytest.approx(math.tanh(3), abs=2e-15)
    assert abs(measured-1) < 0.02


@pytest.mark.parametrize("angle", [30, 45])
@pytest.mark.parametrize("guide", [0, 1])
def test_rotated_harris_2d(estimator, angle, guide):
    measured = measure(estimator, dim=2, angle=angle, guide=guide)
    normal_projection = max(math.cos(math.radians(angle)), math.sin(math.radians(angle)))
    assert measured == pytest.approx(math.tanh(3*normal_projection), abs=2e-15)
    assert abs(measured-1) < 0.03


def test_body_diagonal_harris_3d(estimator):
    measured = measure(estimator, dim=3, guide=1)
    assert measured == pytest.approx(math.tanh(math.sqrt(3)), abs=2e-15)
    assert abs(measured-1) < 0.07


@pytest.mark.parametrize("shear", [60, 120, 180])
def test_rotational_sheet_reconnecting_component(estimator, shear):
    measured = measure(estimator, profile="rotation", shear=shear, guide=1)
    target = math.sin(0.5*math.radians(shear))
    assert measured == pytest.approx(math.sin(0.5*math.radians(shear)*math.tanh(3)),
                                     abs=2e-15)
    assert abs(measured/target-1) < 0.02


@pytest.mark.parametrize("dim", [1, 2, 3])
def test_floor_suppresses_small_grid_noise(estimator, dim):
    assert measure(estimator, profile="noise", dim=dim, amp=0.005,
                   guide=1, floor=0.05) == 0.05


def test_radius_uses_ceiling_on_each_axis(estimator):
    radius = 3.01
    spacing = (0.11, 0.13, 0.17)
    measured = measure(estimator, dim=3, radius=radius, spacing=spacing)
    distance = max(math.ceil(radius/h)*h for h in spacing)/math.sqrt(3)
    assert measured == pytest.approx(math.tanh(distance), abs=2e-15)


@pytest.mark.parametrize("dim", [1, 2, 3])
@pytest.mark.parametrize("nghost", [2, 4])
def test_cached_coefficient_indices_and_constant_fallback(estimator, dim, nghost):
    result = subprocess.run([str(estimator), "cache", str(dim), str(nghost)],
                            capture_output=True, text=True, check=True)
    assert list(map(float, result.stdout.split())) == [0.0, 7.0]
