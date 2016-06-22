//--- Aurora/AuroraLib/RBuffer2D.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/RBuffer2D.hpp"

#include "AuroraLib/Math.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Aurora
{

RBuffer2D::Info RBuffer2D::info() const
{
  const size_t numElements = this->numElements();
  const Real* const data = this->data();

  Info info;

  info.min =  std::numeric_limits<Real>::infinity();
  info.max = -std::numeric_limits<Real>::infinity();
  info.minNorm =  std::numeric_limits<Real>::infinity();
  info.maxNorm = -std::numeric_limits<Real>::infinity();

  for (size_t i = 0; i < numElements; ++i)
  {
    info.min = std::min(info.min, data[i]);
    info.max = std::max(info.max, data[i]);
    info.minNorm = std::min(info.minNorm, powerOf<2>(data[i]));
    info.maxNorm = std::max(info.maxNorm, powerOf<2>(data[i]));
  }

  info.minAbs = std::sqrt(info.minNorm);
  info.maxAbs = std::sqrt(info.maxNorm);

  return info;
}

void RBuffer2D::save(const std::string& filename, const Buffer2DInfo& info,
                     const std::string& description) const
{
  CBuffer2D temp(cols(), rows());
  temp.assign(*this);
  temp.save(filename, info, description);
}

} // namespace Aurora
