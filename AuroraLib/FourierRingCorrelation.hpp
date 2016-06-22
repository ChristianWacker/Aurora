//--- Aurora/AuroraLib/FourierRingCorrelation.hpp ------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_FOURIER_RING_CORRELATION_HPP
#define AURORA_AURORA_LIB_FOURIER_RING_CORRELATION_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <vector>

namespace Aurora
{

/// Class to calculate the Fourier Ring Correlation (FRC) of the absolute
/// values of two ComplexBuffer2D. For the algorithm and a description see
/// [Jensen2010], p. 78f.
class AURORA_API FourierRingCorrelation
{
public:
  /// The Fourier Shell Correlation is defined as [Jensen2010], p.80, eq. 3.5:
  ///   @f[
  ///     FSC(u, v; s) = \frac{\sum_{||s_k|-s|\leq\varepsilon} U(s_k)
  ///       V^*(s_k)} {\left[\left(\sum_{||s_k|-s|\leq\varepsilon}
  ///       \left|U(s_k)\right|^2\right)
  ///       \left(\sum_{||s_k|-s|\leq\varepsilon}\left|V(s_k)\right|^2
  ///       \right)\right]^{1/2}},
  ///   @f]
  /// where @f$s@f$ is the magnitude of the spatial frequency, @f$u@f$ and and
  /// @f$v@f$ are CBuffers2D, @f$U@f$ and @f$V@f$ are their respective
  /// Fourier transformations and @f$2\,\varepsilon@f$ is the thickness of a
  /// ring. In our case the thickness is given by @p binSize.
  /// @param buffer1
  ///  Pointer to the first buffer, i.e. @f$u@f$.
  /// @param buffer2
  ///  Pointer to the second buffer, i.e. @f$v@f$.
  /// @param deltaKX
  ///  Width of a pixel in Fourier space.
  /// @param deltaKY
  ///  Height of a pixel in Fourier space.
  /// @param binSize
  ///  The thickness of a ring @f$2\,\varepsilon@f$.
  /// @remark
  ///  The two buffers must be compatible.
  /// @remark
  ///  For non-square buffer, circles are used. However, depending on the
  ///  source of the buffers, ellipses would probably be more appropriate.
  FourierRingCorrelation(const CBuffer2D& buffer1, const CBuffer2D& buffer2,
                         Real deltaKX, Real deltaKY, Real binSize);

  Real bin(size_t i) const
  {
    assert(i < mNumBins);
    return mResult[i];
  }

  size_t numBins() const
  {
    return mNumBins;
  }

private:
  size_t mNumBins;
  std::vector<Real> mResult;
};

} // namespace Aurora

#endif
