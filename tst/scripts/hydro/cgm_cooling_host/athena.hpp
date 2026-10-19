#ifndef ATHENA_HPP_
#define ATHENA_HPP_
//! \file athena.hpp
//! \brief Host-only stand-in for src/athena.hpp, so src/srcterms/cgm_cooling.hpp builds
//! with a plain C++ compiler (no Kokkos) in the one-zone cooling test.

#include <math.h>

using Real = double;
#define KOKKOS_INLINE_FUNCTION inline

#endif // ATHENA_HPP_
