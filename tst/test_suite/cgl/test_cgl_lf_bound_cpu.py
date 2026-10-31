"""Independent final-stencil LF Jacobian, in temperature coordinates.

Density, staggered B, sound speed, scattering and cap scales are frozen.
VL4 and the complete perpendicular grad-B driver are differentiated. This
does not bound the full nonlinear pressure Jacobian or prove RKL stability.
"""

import json
import os
from pathlib import Path
import shlex
import subprocess

import numpy as np
import pytest


@pytest.mark.parametrize("precision", ["double", "float"])
def test_face_state_extremes(tmp_path, precision):
    # Compile the actual anonymous-TU helper; scattering is disabled in this
    # check, so its unrelated EOS callback can be replaced by zero.
    source = (Path(__file__).resolve().parents[3] / "src/diffusion/cgl_landau_fluid.cpp").read_text()
    structure = source[source.index("struct CGLLFFaceState {"):]
    structure = structure[:structure.index("};")+2]
    start = source.index("bool BuildCGLLFFaceState(")
    function = source[start:source.index("KOKKOS_INLINE_FUNCTION\nvoid CGLLFFlux(", start)]
    unit = tmp_path / "speed.cpp"
    unit.write_text("""
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
using Real = PRECISION;
namespace Kokkos { using std::isfinite; }
struct EOS_Data { Real dfloor=0, pfloor=0, bfloor=0, tfloor=0, nu_coll=0; };
namespace cgl {
Real LimiterCollisionRate(Real, Real, Real, const EOS_Data &, bool) { return 0; }
}
""".replace("PRECISION", precision)+structure+function+"""
int main() {
  const Real minimum=std::numeric_limits<Real>::min();
  const Real maximum=std::numeric_limits<Real>::max();
  const Real subnormal=std::numeric_limits<Real>::denorm_min();
  const Real cases[][3] = {
    {1, 1, 0}, {maximum/4, minimum, 0}, {minimum, maximum/4, 0},
    {minimum, (minimum/subnormal)*(Real(4)/Real(3)), 0},
    {minimum, maximum/4, Real(1e-12)}
  };
  for (const auto &test : cases) {
    EOS_Data eos{}; eos.tfloor=test[2];
    CGLLFFaceState face{};
    if (!BuildCGLLFFaceState(test[1],test[1],test[0],test[0],test[0],test[0],
                            1,0,0,1,0,1,true,1,false,eos,face)) return 1;
    const long double thermal=std::sqrt(static_cast<long double>(test[0]))/
                              std::sqrt(static_cast<long double>(test[1]));
    const long double floor=std::sqrt(static_cast<long double>(test[2]));
    const long double expected=thermal>floor ? thermal : floor;
    const long double relative=std::abs(face.cparallel/expected-1);
    if (!(relative < 8*std::numeric_limits<Real>::epsilon())) {
      std::cerr << "sound speed relative error " << relative << '\\n';
      return 2;
    }
  }
  // Preserve cancellation attenuation and bound CT directions, including norms
  // which overflow Real even though every supplied field component is finite.
  const Real fields[][4] = {
    {Real(.3),Real(.4),0,1}, {1,-2,2,Real(1e-4)},
    {maximum,maximum,maximum,maximum},
    {maximum,-maximum/4,maximum/8,minimum}, {minimum,minimum,0,1}
  };
  for (const auto &field : fields) {
    EOS_Data eos{};
    CGLLFFaceState face{};
    if (!BuildCGLLFFaceState(1,1,1,1,1,1,field[0],field[1],field[2],
                            field[3],0,1,false,1,false,eos,face)) return 3;
    // Scale before normalizing: the true norm can exceed even long double on arm64.
    const long double scale=std::max({std::abs(field[0]),std::abs(field[1]),
                                      std::abs(field[2]),field[3]});
    const long double norm=std::hypot(field[0]/scale,field[1]/scale,field[2]/scale);
    const long double denom=std::max(norm,field[3]/scale);
    const Real actual[3]={face.bhx,face.bhy,face.bhz};
    for (int d=0;d<3;++d) {
      if (!std::isfinite(actual[d]) ||
          std::abs(actual[d]-(field[d]/scale)/denom) >
              8*std::numeric_limits<Real>::epsilon()) return 4;
    }
    if (face.bmag_inv != Real(1)/field[3]) return 5;
  }
}
""")
    executable = tmp_path / "speed"
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    for command in [compiler+["-std=c++17", "-O2", str(unit), "-o", str(executable)],
                    [str(executable)]]:
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout+result.stderr


def test_grad_b_arithmetic_matches_flux(tmp_path):
    # At large anisotropy an ULP of canceled grad-B can dominate the row.
    source = (Path(__file__).resolve().parents[3] / "src/diffusion/cgl_landau_fluid.cpp").read_text()
    start = source.index("Real ScaledMagneticMagnitude(")
    magnitude = source[start:source.index("// Coefficient one-norm", start)]
    start = source.index("Real CGLLFVL4DerivativeNorm(")
    derivative = source[start:source.index("// Positive sums", start)]
    timestep = source[source.index("void CGLLandauFluid::NewTimeStep("):]
    start = timestep.index("const Real delta_b =")
    bound = timestep[start:timestep.index("const Real pressure_ratio =", start)]
    start = source.index("const Real grad_b =", source.index("void CGLLFFlux("))
    flux = source[start:source.index(";", start)+1]
    unit = tmp_path / "grad_b.cpp"
    unit.write_text("""
#include <cmath>
#include <limits>
#include <iostream>
#define KOKKOS_INLINE_FUNCTION
#define SQR(x) ((x)*(x))
using Real=double;
namespace Kokkos { using std::isfinite; }
struct Face { Real bhx,bhy,bhz,bhdir; };
"""+magnitude+derivative+"""
Real field(Real x,Real y) { return 1+.24*(x-y); }
Real bcc(int i,int j,int d) {
  return d==0 ? .5*field(i-.5,j)+.5*field(i+.5,j) :
                .5*field(i,j-.5)+.5*field(i,j+.5);
}
Real magnitude(int i,int j) {
  return ScaledMagneticMagnitude(bcc(i,j,0),bcc(i,j,1),0);
}
int main() {
  int failed=0;
  for (int scenario=0;scenario<2;++scenario) {
    // First case: divergence-free affine CT field; second isolates 3D sum order.
    const int ndim=scenario==0 ? 2 : 3, dir=scenario==0 ? 0 : 2;
    const int m=0,il=0,jl=0,kl=0,ir=dir==0,jr=0,kr=dir==2;
    const Real bar=.5*magnitude(0,0)+.5*magnitude(1,0);
    const Real bhx=field(.5,0)/bar, bhy=(.5*bcc(0,0,1)+.5*bcc(1,0,1))/bar;
    const Face face=scenario==0 ? Face{bhx,bhy,0,bhx} : Face{1,1,1,1};
    const Real dx[3]={scenario==0 ? .09 : 1, scenario==0 ? .09 : 1,4};
    const Real bh[3]={face.bhx,face.bhy,face.bhz}, bv[3]={bh[0],bh[1],bh[2]};
    const Real minimum=std::numeric_limits<Real>::min();
    Real grad_parallel=0,grad_perp=0,derivative_parallel[3]={},derivative_perp[3]={};
    auto bmag=[&](int,int k,int j,int i) {
      return scenario==0 ? magnitude(i,j) : 3e16+1e16*i-1e16*j+4*k;
    };
    auto tpar=[](int,int,int,int) { return Real(1); };
    auto tperp=tpar;
    const Real selected=[&]() {
"""+bound+"""
      return grad_b;
    }();
    Real gradient[3]={};
    for (int t=0;t<ndim;++t) {
      const int di=t==0,dj=t==1,dk=t==2;
      gradient[t]=t==dir ? (bmag(m,kr,jr,ir)-bmag(m,kl,jl,il))/dx[t] :
          .25*(bmag(m,kr+dk,jr+dj,ir+di)-bmag(m,kr-dk,jr-dj,ir-di)+
               bmag(m,kl+dk,jl+dj,il+di)-bmag(m,kl-dk,jl-dj,il-di))/dx[t];
    }
    const Real gb_x=gradient[0],gb_y=gradient[1],gb_z=gradient[2];
"""+flux+"""
    if (selected!=grad_b || grad_b==0) {
      std::cerr << "scenario " <<scenario<< ": bound " <<selected
                << ", flux " <<grad_b<<'\\n';
      failed |= 1<<scenario;
    }
  }
  return failed;
}
""")
    executable = tmp_path / "grad_b"
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    for command in [compiler+["-std=c++17", "-O2", "-ffp-contract=off", str(unit), "-o", str(executable)],
                    [str(executable)]]:
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout+result.stderr


@pytest.mark.parametrize("normal_field", [.1, 10.])
def test_cold_grad_b_envelope_matches_flux(tmp_path, normal_field):
    root = Path(__file__).resolve().parents[3]
    source = (root / "src/diffusion/cgl_landau_fluid.cpp").read_text()
    physics = (root / "src/eos/cgl_physics.hpp").read_text()
    eos = (root / "src/eos/ideal_c2p_mhd.hpp").read_text()

    def function(text, signature):
        start = text.index("\n"+signature)+1
        end, depth = text.index("{", start)+1, 1
        while depth:
            depth += (text[end] == "{") - (text[end] == "}")
            end += 1
        return text[start:end]+"\n"

    timestep = source[source.index("void CGLLandauFluid::NewTimeStep("):]
    gradient = timestep[timestep.index("const Real delta_b ="):
                        timestep.index("const Real pressure_ratio =")]
    cold = timestep[timestep.index("const Real log_bar ="):
                    timestep.index("log_perp + CGLLFLogAbs(bi) - log_bar;")+
                    len("log_perp + CGLLFLogAbs(bi) - log_bar;")]
    unit = tmp_path / "cold_grad_b.cpp"
    unit.write_text(("""
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#define SQR(x) ((x)*(x))
using Real=double;
namespace Kokkos {
using std::isfinite; using std::isinf; using std::isnan;
using std::log2; using std::exp2;
}
namespace cgl {
"""+physics[physics.index("constexpr Real kSqrtTwoOverPi"):
             physics.index("KOKKOS_INLINE_FUNCTION")]+
        function(physics, "Real PerpendicularHeatFluxRatio(")+"}\n"+
        "".join(function(source, signature) for signature in [
            "Real ScaledMagneticMagnitude(", "Real CGLLFVL4DerivativeNorm(",
            "Real CGLLFLogAbs(", "Real CGLLFLogAdd(", "Real CGLLFLogScale(",
            "Real CGLLFLogDiffusivity("])+
        function(eos, "Real CGLConservedAnisotropy(")+
        function(eos, "void CGLRecoverPressuresFromInternalEnergyAndAnisotropy(")+"""
int main() {
  const Real h=1e-160,cp=1e-310,kpar=1e307,nu=0,rho=1;
  Real z=std::sqrt(.98),za=z;
  for (int n=0;n<6;++n) za=std::nextafter(za,INFINITY);
  // Begin with a canceled grad-B projection, then also cover a CT normal
  // component larger than the cell-centered magnetic magnitudes.
  const Real a=ScaledMagneticMagnitude(.1,.1,za);
  const Real b=ScaledMagneticMagnitude(.1,.1,z);
  const Real c=ScaledMagneticMagnitude(.1,.1,std::sqrt(16-.02));
  const Real bbar=b,bi=b,beta=bi/bbar;
  const Real normal_field=NORMAL_FIELD;
  const Real direction=std::fmax(bbar,std::sqrt(normal_field*normal_field+.01));
  const Real bn=normal_field/direction,bt=.1/direction;
  Real pp,pt;
  const Real magnetic_energy=.5*(.02+z*z);
  const Real energy=1+magnetic_energy;
  const Real anisotropy=CGLConservedAnisotropy(1,1e-160,1,b);
  CGLRecoverPressuresFromInternalEnergyAndAnisotropy(
      1,energy-magnetic_energy,anisotropy,b,pp,pt);
  if (!(pp>0 && pt>0)) return 1;
  struct Face { Real bhx,bhy,bhz,bhdir,ppar,pperp,rho; };
  const Face face{bn,bt,0,bn,pp,pt,1};
  const Real dx[3]={h,h,1},bh[3]={bn,bt,0};
  Real bv[3]={normal_field,.1,0};
  const int m=0,kl=0,jl=0,il=0,kr=0,jr=0,ir=1,ndim=2,dir=0;
  const Real minimum=std::numeric_limits<Real>::min();
  Real grad_parallel=0,grad_perp=0,derivative_parallel[3]={},derivative_perp[3]={};
  auto bmag=[&](int,int,int j,int i) {
    if (i==0 && std::abs(j)==1) return c;
    return i==1 && j==1 ? a : b;
  };
  auto tpar=[&](int,int,int j,int i) { return j==0 && (i==0 || i==1) ? pp : Real(1); };
  auto tperp=[&](int,int,int,int) { return pt; };
"""+gradient+"""
  auto row=[&](Real represented_gradient) {
    const Real grad_b=represented_gradient;
    Real log_face_parallel,log_face_perp;
"""+cold+"""
    return std::exp2(log_face_perp);
  };
  const Real bound=row(grad_b);
  const Real qratio=cgl::PerpendicularHeatFluxRatio(cp,1,pp,pt,1/bbar,kpar,0,0,grad_b);
  const long double ratio=static_cast<long double>(pt)/pp;
  const long double chi=static_cast<long double>(cgl::kSqrtTwoOverPi)*cp/kpar;
  // The parallel-to-perpendicular block alone is a lower bound on the full row.
  const long double reverse=chi*bn*ratio*ratio*std::abs(grad_b)/bbar/h;
  const long double capped=reverse/((1+std::abs(qratio))*(1+std::abs(qratio)));
  std::cout<<std::setprecision(17)<<"cold row "<<bound
           <<", actual capped reverse row "<<capped<<", qratio "<<qratio<<'\\n';
  if (!(std::isfinite(bound) && bound>=capped*(1-3e-12))) {
    std::cerr<<"cold row "<<bound<<", actual capped reverse row "<<capped<<'\\n';
    return 2;
  }
  // Preserve the absolute-sum estimate if the represented projection vanishes;
  // invalid projections must never silently select a finite row.
  if (!(row(0)>0 && std::isfinite(row(0)))) return 3;
  if (!std::isinf(row(INFINITY))) return 4;
  if (!std::isnan(row(std::numeric_limits<Real>::quiet_NaN()))) return 5;
  // Once the direction is bounded, scaling the entire face vector changes
  // neither its direction nor the cold row; Bbar still controls focusing.
  if (normal_field>bbar) {
    bv[0]*=2; bv[1]*=2;
    if (std::abs(row(grad_b)/bound-1)>3e-12) return 6;
  }
}
""").replace("NORMAL_FIELD", repr(normal_field)))
    executable = tmp_path / "cold_grad_b"
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    for command in [compiler+["-std=c++17", "-O2", "-ffp-contract=off", str(unit), "-o", str(executable)],
                    [str(executable)]]:
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout+result.stderr


class _FrozenStencil:
    def __init__(self, rho, face_b, state, spacing, nu=0.0, kpar=2*np.pi,
                 local=False, cap=True):
        self.shape, self.n = rho.shape, rho.size
        self.ndim, self.h = rho.ndim, np.asarray(spacing)
        self.rho, self.state = rho.ravel(), np.asarray(state).reshape(2, -1)
        self.coords = list(np.ndindex(self.shape))
        bcc = np.array(face_b, copy=True)
        for d in range(self.ndim):
            bcc[d] = (face_b[d]+np.roll(face_b[d], 1, axis=d))/2
        bcc = bcc.reshape(3, -1)
        self.bmag = np.linalg.norm(bcc, axis=0)
        self.faces = []
        for d in range(self.ndim):
            for c in self.coords:
                l, r = self.index(c), self.index(self.move(c, d, 1))
                bar = (self.bmag[l]+self.bmag[r])/2
                if bar == 0:
                    continue
                bv = (bcc[:, l]+bcc[:, r])/2
                bv[d] = face_b[(d,)+c]
                bh = bv/max(bar, np.linalg.norm(bv))
                rf = (self.rho[l]+self.rho[r])/2
                mean = np.zeros(self.n)
                mean[[l, r]] = self.rho[[l, r]]/(2*rf)
                tp, tt = self.state @ mean
                speed = np.sqrt(tp) if local else 1.0
                chi = np.array([
                    8*speed**2/(np.sqrt(8*np.pi)*speed*kpar+(3*np.pi-8)*nu),
                    speed**2/(np.sqrt(np.pi/2)*speed*kpar+nu),
                ])
                gb = sum(bh[t]*self.gradient(self.bmag, c, d, t, False)[0]
                         for t in range(self.ndim))/bar
                qmax = speed*rf*np.array([np.sqrt(8/np.pi)*tp, np.sqrt(2/np.pi)*tt])
                if not cap:
                    qmax[:] = np.inf
                self.faces.append((c, d, l, r, bh, rf, mean, chi, gb, qmax, bar))

    def index(self, c):
        return np.ravel_multi_index(c, self.shape)

    def move(self, c, d, offset):
        result = list(c)
        result[d] = (result[d]+offset) % self.shape[d]
        return tuple(result)

    def gradient(self, state, c, normal, direction, limited=True):
        """Value, exact branch derivative, and generalized row envelope."""
        right = self.move(c, normal, 1)
        l, r = self.index(c), self.index(right)
        row = np.zeros(self.n)
        if direction == normal:
            row[[l, r]] = [-1, 1]
            return ((state[r]-state[l])/self.h[direction], row/self.h[direction],
                    2/self.h[direction])
        rp, rm, lp, lm = (self.index(self.move(cell, direction, step))
                          for cell, step in [(right, 1), (right, -1), (c, 1), (c, -1)])
        slopes = np.array([state[rp]-state[r], state[r]-state[rm],
                           state[lp]-state[l], state[l]-state[lm]])
        partial = np.zeros(4)
        value = 0.0
        if not limited:
            value, partial[:] = np.mean(slopes), 0.25
        elif np.all(slopes > 0) or np.all(slopes < 0):
            value = 4/np.sum(1/slopes)
            partial = value**2/(4*slopes**2)
        for index, coefficient in [(rp, partial[0]), (r, partial[1]-partial[0]),
                                   (rm, -partial[1]), (lp, partial[2]),
                                   (l, partial[3]-partial[2]), (lm, -partial[3])]:
            row[index] += coefficient
        envelope = 2*(max(partial[:2])+max(partial[2:]))
        if limited and np.any(slopes == 0) and (np.all(slopes >= 0) or np.all(slopes <= 0)):
            envelope = 8.0
        return value/self.h[direction], row/self.h[direction], envelope/self.h[direction]

    def evaluate(self, state=None, matrix=False):
        state = self.state if state is None else np.asarray(state).reshape(2, -1)
        result = np.zeros((2, self.n))
        jacobian = np.zeros((2*self.n, 2*self.n))
        bound = np.zeros((2, self.n))
        for c, d, l, r, bh, rf, mean, chi, gb, qmax, bar in self.faces:
            grad = np.zeros(2)
            derivative = np.zeros((2, self.n))
            envelope = np.zeros(2)
            for channel in range(2):
                for t in range(self.ndim):
                    value, row, norm = self.gradient(state[channel], c, d, t)
                    grad[channel] += bh[t]*value
                    derivative[channel] += bh[t]*row
                    envelope[channel] += abs(bh[t])*norm
            tp, tt = state @ mean
            ratio = tt/tp
            grad[1] -= tt*(1-ratio)*gb
            raw = -rf*chi*grad
            cap = 1/(1+abs(raw)/qmax)
            flux = raw*cap*bh[d]
            if matrix:
                derivative[1] -= (1-2*ratio)*gb*mean
                reverse = -ratio**2*gb*mean
                flux_rows = np.zeros((2, 2*self.n))
                flux_rows[0, :self.n] = derivative[0]
                flux_rows[1, :self.n] = reverse
                flux_rows[1, self.n:] = derivative[1]
                flux_rows *= (-rf*chi*cap**2*bh[d])[:, None]
                envelope[1] += (abs(1-2*ratio)+ratio**2)*abs(gb)
                envelope *= chi*abs(bh[d])
            for i, sign in [(l, -1), (r, 1)]:
                factor = sign/(self.rho[i]*self.h[d])
                beta = self.bmag[i]/bar
                result[:, i] += factor*np.array([flux[0]+2*(1-beta)*flux[1], beta*flux[1]])
                if matrix:
                    jacobian[i] += factor*(flux_rows[0]+2*(1-beta)*flux_rows[1])
                    jacobian[self.n+i] += factor*beta*flux_rows[1]
                    bound[:, i] += rf*abs(factor)*np.array([
                        envelope[0]+2*abs(1-beta)*envelope[1], beta*envelope[1]])
        return (jacobian, bound.ravel()) if matrix else result.ravel()

    def matrix_and_bound(self):
        return self.evaluate(matrix=True)


# The retained trial 36 disproves the old secant/coefficient-five envelope.
_COUNTEREXAMPLE = np.array([
    -.2931815757973079, 1.3265700611876208, 1.090907140691816, -1.8315350368603842,
    -.43873315311167316, 1.1063250396732305, -.9842393273663874, -.3537786596348982,
    -.15619196817321634, -.7757273661906982, .7843296665164378, 1.396450244892953,
    -.34746366904753334, -1.0645902368437121, .6959146227470747, -.8175773742732974,
    .7672748172933098, .29023923882679087, .5527340239744399, 1.068899314299884,
    -1.3493815525616166, -1.3410664127000758, .32372521336423055, -.6832774689277039,
    .9607161230179032, -1.1429005853858107, -.40945921621049, .5240165403444472,
    -.8082724613006674, 1.7766746906767819, -2.0149880625585985, -.3794956345669529,
    -.37079306846264376, -.10182373194345405, -.6334832636400648, .841314905834914,
])


def _case(name):
    shape = (3, 3, 3) if name == "3d" else ((6, 6) if name == "counterexample" else (4, 4))
    coordinates = np.meshgrid(*[(np.arange(n)+.5)/n for n in shape], indexing="ij")
    rho = np.ones(shape)
    face_b = np.zeros((3, *shape))
    face_b[0], face_b[1] = np.cos(np.pi/6), np.sin(np.pi/6)
    if name == "3d":
        face_b[2] = .4
    if name == "checkerboard":
        alternating = (-1.)**np.indices(shape).sum(axis=0)
        face_b[0], face_b[1] = .1+.9*alternating, .1-.9*alternating
        divergence = sum(n*(face_b[d]-np.roll(face_b[d], 1, axis=d))
                         for d, n in enumerate(shape))
        np.testing.assert_array_equal(divergence, 0)
    if name == "density":
        rho = np.where(coordinates[0] < .5, 1., 200.)
    if name == "gradb":
        face_b[0] = .03
        face_b[1] = np.tanh(np.sin(2*np.pi*coordinates[0])/.5)
        # Break the four-cell magnitude degeneracy while keeping div B = 0.
        face_b[1] += .27
    phase = sum((d+1.137)*v for d, v in enumerate(coordinates))
    state = np.array([1+.07*np.sin(2*np.pi*phase+.13),
                      .7+.04*np.cos(2*np.pi*phase+.37)])
    if name == "counterexample":
        face_b[0] = face_b[1] = np.sqrt(.5)
        state[0] = 1+.01*_COUNTEREXAMPLE.reshape(shape)
        state[1] = 1
    spacing = [1/n for n in shape]
    if name == "unequal":
        spacing = [.13, .31]
    return _FrozenStencil(rho, face_b, state, spacing, nu=.8 if name == "gradb" else 0.,
                          local=name == "gradb", cap=name != "counterexample")


def _rkl_amplification(z, stages):
    """Legendre polynomial, independently of the production stage recurrence."""
    denominator = stages*stages+stages-2
    x = 1+4*z/denominator
    before, current = np.ones_like(x), x
    for degree in range(2, stages+1):
        before, current = current, ((2*degree-1)*x*current-(degree-1)*before)/degree
    weight = denominator/(2*stages*(stages+1))
    return 1-weight+weight*current


@pytest.mark.parametrize("name", ["oblique", "unequal", "density", "checkerboard",
                                  "gradb", "3d", "counterexample"])
def test_frozen_final_stencil(name, tmp_path):
    grid = _case(name)
    matrix, envelope = grid.matrix_and_bound()
    finite_difference = np.empty_like(matrix)
    for j in range(2*grid.n):
        perturbation = np.zeros(2*grid.n)
        perturbation[j] = 2e-7
        finite_difference[:, j] = (grid.evaluate(grid.state.ravel()+perturbation)
                                    -grid.evaluate(grid.state.ravel()-perturbation))/(4e-7)
    scale = max(abs(matrix).sum(axis=1))
    assert max(abs(matrix-finite_difference).sum(axis=1)) < 2e-6*scale
    rows = abs(matrix).sum(axis=1)
    assert np.all(rows <= envelope+1e-12*max(envelope))
    energy_weights = np.r_[.5*grid.rho, grid.rho]
    np.testing.assert_allclose(energy_weights @ matrix, 0, rtol=0,
                               atol=1e-13*max(energy_weights)*scale)
    eigenvalues = np.linalg.eigvals(matrix)
    assert max(abs(eigenvalues)) <= max(envelope)*(1+1e-12)
    if name == "gradb":
        assert np.max(abs(matrix[grid.n:, :grid.n])) > 100*np.finfo(float).eps*scale
    if name == "counterexample":
        # Same retained state at a larger amplitude avoids finite-difference
        # cancellation; VL4 partials are invariant under positive rescaling.
        assert max(abs(eigenvalues)) > 1.25*54.85850144718421
        assert max(eigenvalues.real) > .1

    safety = .3
    explicit_dt = safety*2/max(envelope)
    sweep_dt = 12*explicit_dt
    stages = 3
    while (stages*stages+stages-2)/4 < sweep_dt/explicit_dt:
        stages += 2
    amplification = _rkl_amplification(sweep_dt*eigenvalues, stages)
    exact = np.exp(sweep_dt*eigenvalues)
    # Positive/complex eigenvalues may represent actual frozen-operator growth.
    # Check time refinement toward the explicit reference, not |R| <= 1.
    errors = []
    for subdivisions in [1, 2, 4]:
        result = _rkl_amplification(sweep_dt*eigenvalues/subdivisions, stages)**subdivisions
        errors.append(float(max(abs(result-exact))))
    assert errors[2] < errors[1] < errors[0]
    report = dict(name=name, row_norm=float(max(rows)), envelope=float(max(envelope)),
                  explicit_dt=explicit_dt, sweep_dt=sweep_dt, stages=stages,
                  eigenvalues=np.column_stack((eigenvalues.real, eigenvalues.imag)).tolist(),
                  amplification_abs=abs(amplification).tolist(),
                  exact_amplification_abs=abs(exact).tolist(), refinement_errors=errors)
    (tmp_path / "frozen-stencil.json").write_text(json.dumps(report, indent=2)+"\n")
