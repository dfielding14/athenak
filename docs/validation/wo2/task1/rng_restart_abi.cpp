#include <cstddef>
#include "utils/random.hpp"
#include "srcterms/turb_driver.hpp"
static_assert(sizeof(RNG_State) == 296);
static_assert(offsetof(RNG_State, idum) == 0);
static_assert(offsetof(RNG_State, idum2) == 8);
static_assert(offsetof(RNG_State, iy) == 16);
static_assert(offsetof(RNG_State, iv) == 24);
static_assert(offsetof(RNG_State, iset) == 280);
static_assert(sizeof(RNG_State::iset) == 4);
static_assert(offsetof(RNG_State, gset) == 288);
static_assert(sizeof(RNG_State::gset) == 8);
static_assert(sizeof(TurbulenceRestartMetadata) == 248);
