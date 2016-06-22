//--- Aurora/AuroraLib/Laplace.cpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Laplace.hpp"

namespace Aurora
{

#if (AURORA_COMPILER == AURORA_COMPILER_MSVC)
const Real Laplace<ThreePoint>::factor[] = {-4.0, 1.0};
const Real Laplace<FivePoint> ::factor[] = {-5.0, 4.0 / 3.0, -1.0 / 12.0};
const Real Laplace<SevenPoint>::factor[]
  = {-49.0 / 9.0, 1.5, -0.15, 1.0 / 90.0};
const Real Laplace<NinePoint> ::factor[]
  = {-205.0 / 36.0, 1.6, -0.2, 8.0 / 315.0, -1.0 / 560.0};
#else
constexpr Real Laplace<ThreePoint>::factor[];
constexpr Real Laplace<FivePoint> ::factor[];
constexpr Real Laplace<SevenPoint>::factor[];
constexpr Real Laplace<NinePoint> ::factor[];
#endif

} // namespace Aurora
