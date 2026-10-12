//========================================================================================
// AthenaK astrophysical plasma code
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \brief Independent-kernel regression for passive wall invariant decoding.

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_passive.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"

namespace {
void Require(const char *name, bool condition) {
  if (!condition) {
    std::cout << "Passive wall decode regression failed: " << name << std::endl;
    std::exit(EXIT_FAILURE);
  }
}
}  // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  Require("double precision", sizeof(Real) == sizeof(double));
  auto *pack = pmy_mesh_->pmb_pack;
  auto *pmhd = pack->pmhd;
  auto *peos = pmhd->peos;
  const auto eos = peos->eos_data;
  Require("passive no-backup model", eos.passive && !eos.backup_lim);
  auto cons = pmhd->u0;
  auto prim = pmhd->w0;
  auto bcc = pmhd->bcc0;
  auto b = pmhd->b0;
  auto copy_to_host = [](const DvceArray5D<Real> &view) {
    auto host = Kokkos::create_mirror(view);
    Kokkos::deep_copy(host, view);
    return host;
  };
  const auto indcs = pmy_mesh_->mb_indcs;
  const int n1 = cons.extent_int(4), n2 = cons.extent_int(3), n3 = cons.extent_int(2);
  const int nmb = pack->nmb_thispack;
  // Exact wall-encoded state from strict HIP preflight, cycle 69, t=0.15569028418406383.
  // Previously the actual grid decoder gave pperp-ppar+B^2=-6.661338147750939e-16.
  const Real bx = -0.087016078880545439, by = 0.0375499263781175;
  const Real bz = 0.93842460882593404;
  const Real captured_rho = 1.0507151506792887;
  const Real captured_j = 1.5868588473676859, captured_a = 0.11669182471092962;
  Kokkos::deep_copy(b.x1f, bx);
  Kokkos::deep_copy(b.x2f, by);
  Kokkos::deep_copy(b.x3f, bz);
  Kokkos::deep_copy(cons, 0.0);
  par_for("passive_wall_fixture", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    cons(m,IDN,k,j,i) = captured_rho;
    cons(m,IEN,k,j,i) = captured_j;
    cons(m,IAN,k,j,i) = captured_a;
  });
  // This is the production C2P kernel, independently compiled from the encoder below.
  peos->ConsToPrim(cons, b, prim, bcc, false, 0,n1-1, 0,n2-1, 0,n3-1);
  auto check = [&](const char *label, const bool must_be_admissible = true) {
    int bad = 0;
    Kokkos::parallel_reduce("passive_wall_independent_check", nmb*n1*n2*n3,
    KOKKOS_LAMBDA(int index, int &nbad) {
      const int i=index%n1, j=(index/n1)%n2, k=(index/(n1*n2))%n3;
      const int m=index/(n1*n2*n3);
      const Real bsqr=SQR(bcc(m,IBX,k,j,i))+SQR(bcc(m,IBY,k,j,i))+SQR(bcc(m,IBZ,k,j,i));
      if (cgl::HardBoundViolated(prim(m,IPP,k,j,i)-prim(m,IPR,k,j,i), bsqr, eos, false)) ++nbad;
    }, Kokkos::Sum<int>(bad));
    auto host = copy_to_host(prim);
    const Real pp=host(0,IPR,indcs.ks,indcs.js,indcs.is);
    const Real pt=host(0,IPP,indcs.ks,indcs.js,indcs.is);
    std::cout << std::setprecision(17) << label << ": ppar=" << pp << " pperp=" << pt
              << " hard_bad=" << bad << std::endl;
    if (must_be_admissible) Require(label, bad == 0);
  };
  // The old payload may still be outside under the new canonical arithmetic.
  // Re-encode it through the actual wall operator; admissibility must then be
  // independent of the encoder's kernel context, as at the real split checkpoint.
  check("captured old encoding before reprojection", false);
  peos->Collisions(prim, bcc, cons, 0.0, CGLCollisionMode::walls_only,
                  0,n1-1, 0,n2-1, 0,n3-1);
  peos->ConsToPrim(cons, b, prim, bcc, false, 0,n1-1, 0,n2-1, 0,n3-1);
  check("captured wall survives canonical C2P");

  // Vary rho/U at the same magnetic field and test both wall-encoded and
  // genuinely outside-wall states through the production collision operator.
  par_for("passive_wall_encoding_kernel", DevExeSpace(), 0,nmb-1, 0,n3-1, 0,n2-1, 0,n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const int n=i+n1*(j+n2*(k+n3*m));
    const Real rho = 0.85 + 0.3*static_cast<Real>((37*n)%1009)/1009.0;
    const Real u = 6.5 + 2.0*static_cast<Real>((71*n)%1013)/1013.0;
    const Real bsqr = SQR(bx)+SQR(by)+SQR(bz);
    Real pp, pt;
    cgl::PassivePair q;
    if (n%2 == 0) {
      q = cgl::PassiveWallEncode(rho, u, -bsqr, bsqr, eos, false, pp, pt);
    } else {
      const Real delta = -1.1*bsqr;
      pp = TWO_3RDS*u - TWO_3RDS*delta;
      pt = TWO_3RDS*u + ONE_3RD*delta;
      q = cgl::PassiveEncode(rho, pp, pt, sqrt(bsqr));
    }
    cons(m,IDN,k,j,i)=rho;
    cons(m,IEN,k,j,i)=q.j;
    cons(m,IAN,k,j,i)=q.a;
  });
  peos->ConsToPrim(cons, b, prim, bcc, false, 0,n1-1, 0,n2-1, 0,n3-1);
  auto before = copy_to_host(cons);
  auto before_w = copy_to_host(prim);
  peos->Collisions(prim, bcc, cons, 0.0, CGLCollisionMode::walls_only,
                  0,n1-1, 0,n2-1, 0,n3-1);
  auto local = copy_to_host(prim);
  peos->ConsToPrim(cons, b, prim, bcc, false, 0,n1-1, 0,n2-1, 0,n3-1);
  check("encoded walls survive independent canonical C2P");
  auto after = copy_to_host(cons);
  auto after_w = copy_to_host(prim);
  for (int m=0; m<nmb; ++m) for (int k=0; k<n3; ++k) {
    for (int j=0; j<n2; ++j) for (int i=0; i<n1; ++i) {
      for (int n=IDN; n<=IM3; ++n) {
        Require("wall preserves flow bits", before(m,n,k,j,i)==after(m,n,k,j,i));
      }
      for (int n : {IPR, IPP}) {
        Require("local and canonical pressure agree", local(m,n,k,j,i)==after_w(m,n,k,j,i));
      }
      const Real u_before=0.5*before_w(m,IPR,k,j,i)+before_w(m,IPP,k,j,i);
      const Real u_after=0.5*after_w(m,IPR,k,j,i)+after_w(m,IPP,k,j,i);
      Require("wall preserves physical thermal energy", std::abs(u_after-u_before)
          <= 64*std::numeric_limits<Real>::epsilon()*u_before);
    }
  }
  // Copy only persistent conserved data, as a restart would, and repeat the
  // separately compiled decoder and collision map. Neither may drift the bits.
  Kokkos::deep_copy(cons, after);
  peos->ConsToPrim(cons, b, prim, bcc, false, 0,n1-1, 0,n2-1, 0,n3-1);
  auto restart_w = copy_to_host(prim);
  Require("persistent reload primitive bits", std::memcmp(after_w.data(), restart_w.data(),
          after_w.size()*sizeof(Real))==0);
  peos->Collisions(prim, bcc, cons, 0.0, CGLCollisionMode::walls_only,
                  0,n1-1, 0,n2-1, 0,n3-1);
  auto repeated = copy_to_host(cons);
  Require("wall fixed point", std::memcmp(after.data(), repeated.data(), after.size()*sizeof(Real))==0);
  std::cout << "Passive production wall decode checks passed" << std::endl;
}
