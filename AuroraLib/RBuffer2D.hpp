//--- Aurora/AuroraLib/RBuffer2D.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_R_BUFFER_2D_HPP
#define AURORA_AURORA_LIB_R_BUFFER_2D_HPP

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/MarginBuffer2D.hpp"

namespace Aurora
{

/// This class represents a two dimensional array of real values.
class AURORA_API RBuffer2D : public MarginBuffer2D<RBuffer2D, Real, 0>
{
private:
  typedef MarginBuffer2D<RBuffer2D, Real, 0> Base;

public:
  struct Info
  {
    Real min;
    Real max;
    Real minNorm;
    Real maxNorm;
    Real minAbs;
    Real maxAbs;
  };

  RBuffer2D() = default;

  // Use the constructors from the base class
  using MarginBuffer2D<RBuffer2D, Real, 0>::MarginBuffer2D;

  RBuffer2D(RBuffer2D&&) noexcept = default;
  RBuffer2D& operator=(RBuffer2D&&) noexcept = default;

  /// delete the copy constructor
  RBuffer2D(const RBuffer2D&) = delete;

  /// delete the copy assignment operator
  RBuffer2D& operator=(const RBuffer2D&) = delete;

  Info info() const;

  void save(const std::string& filename, const Buffer2DInfo& info,
            const std::string& description) const;
};

} // namespace Aurora

#endif
