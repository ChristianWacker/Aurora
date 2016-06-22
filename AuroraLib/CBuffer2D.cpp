//--- Aurora/AuroraLib/CBuffer2D.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CBuffer2D.hpp"

namespace Aurora
{

template<>
AURORA_API typename CMarginBuffer2D<0>::Loaders& CMarginBuffer2D<0>::loaders()
{
  static Loaders loaders;
  return loaders;
}

template<>
AURORA_API typename CMarginBuffer2D<0>::Savers& CMarginBuffer2D<0>::savers()
{
  static Savers savers;
  return savers;
}

} // namespace Aurora
