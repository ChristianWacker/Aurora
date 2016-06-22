//--- Aurora/AuroraLib/CtfIncoherent.cpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CtfIncoherent.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Physics.hpp"

namespace Aurora
{

CtfIncoherent::CtfIncoherent(Real energy) :
  mCoherent(energy), mCC(0.0), mHighTensionRipple(0.0), mCurrentInstability(0.0),
  mEnergySpread(0.0), mSigmaBeta(0.0)
{ }

CBuffer2D CtfIncoherent::buffer(const Buffer2DInfoBase& info,
                                Real factor) const
{
  const int64_t cols = info.cols();
  const int64_t rows = info.rows();

  CBuffer2D result(cols, rows);

  const Real deltaThetaX = info.deltaKX() / mCoherent.wavenumber();
  const Real deltaThetaY = info.deltaKY() / mCoherent.wavenumber();

  #pragma omp parallel for
  for (int64_t j = 0; j < rows; ++j)
  {
    // y coordinate in Fourier space
    const Real thetaY = deltaThetaY * (j > rows / 2 ? j - rows : j);

    for (int64_t i = 0; i < cols; ++i)
    {
      // x coordinate in Fourier space
      const Real thetaX = deltaThetaX * (i > cols / 2 ? i - cols : i);
      const Real theta = std::hypot(thetaX, thetaY);
      result.pixel(i, j) = factor * mCoherent.value(thetaX, thetaY) *
                           temporalEnvelope(theta) * spatialEnvelope(theta);
    }
  }

  return result;
}

Real CtfIncoherent::temporalEnvelope(Real theta) const
{
  return std::exp(R(-0.125) * powerOf<2>(focusSpread() * mCoherent.wavenumber() *
                  powerOf<2>(theta)));
}

Real CtfIncoherent::spatialEnvelope(Real theta) const
{
  return std::exp(R(-0.5) * powerOf<2>(mSigmaBeta * mCoherent.wavenumber() *
                  (mCoherent.c(1) + mCoherent.c(3) * powerOf<2>(theta)) * theta));
}

} // namespace Aurora
