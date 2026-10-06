#include "athena.hpp"
#include "eos/eos.hpp"
#include "eos/cgl_passive.hpp"
static_assert(sizeof(Real) == sizeof(double), "This proof is for final double builds");
#if defined(KOKKOS_ENABLE_HIP)
#define PROOF_EXPORT __host__ __device__ __attribute__((noinline, used))
#else
#define PROOF_EXPORT __attribute__((noinline, used))
#endif
extern "C" PROOF_EXPORT
void wo2_passive_specific(Real rho, Real ppar, Real pperp, Real bmag, Real *out) {
  const auto q = cgl::PassiveSpecificInvariants(rho, ppar, pperp, bmag);
  out[0] = q.j;
  out[1] = q.a;
}
extern "C" PROOF_EXPORT
void wo2_passive_encode(Real rho, Real ppar, Real pperp, Real bmag, Real *out) {
  const auto q = cgl::PassiveEncode(rho, ppar, pperp, bmag);
  out[0] = q.j;
  out[1] = q.a;
}
#if defined(KOKKOS_ENABLE_HIP)
extern "C" __global__
void wo2_passive_kernel(const Real *values, Real *out) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  wo2_passive_specific(values[4*i], values[4*i+1], values[4*i+2],
                       values[4*i+3], out+4*i);
  wo2_passive_encode(values[4*i], values[4*i+1], values[4*i+2],
                     values[4*i+3], out+4*i+2);
}
#endif
