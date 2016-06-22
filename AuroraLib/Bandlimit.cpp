//--- Aurora/AuroraLib/Bandlimit.cpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Bandlimit.hpp"

#include "AuroraLib/MultisliceParameters.hpp"

namespace Aurora
{

std::unique_ptr<CBandlimit> CBandlimit::create(Bandlimit bandlimit,
                                               const Buffer2DInfo& info,
                                               const CBuffer2D& bufferIn,
                                               CBuffer2D& bufferOut)
{
  switch (bandlimit)
  {
  case Bandlimit::Disabled:
    return std::make_unique<CBandlimitDisabled>(info, bufferIn, bufferOut);
  case Bandlimit::Sharp:
    return std::make_unique<CBandlimitSharp>(info, bufferIn, bufferOut);
  case Bandlimit::Smooth:
    return std::make_unique<CBandlimitSmooth>(info, bufferIn, bufferOut);
  default:
    AURORA_UNREACHABLE;
  }
}

RBuffer2D CBandlimitSharp::mask(const Buffer2DInfo& info)
{
  const int64_t cols = info.cols();
  const int64_t rows = info.rows();
  const Real deltaKX = info.deltaKX();
  const Real deltaKY = info.deltaKY();
  const Real kXNyquistSqr = powerOf<2>(info.kXNyquist());
  const Real kYNyquistSqr = powerOf<2>(info.kYNyquist());
  const Real numPixels = info.numPixels<Real>();

  RBuffer2D result(cols, rows);

  #pragma omp parallel for
  for (int64_t j = 0; j < rows; ++j)
  {
    // coordinates in reciprocal space
    const Real kY = deltaKY * (j > rows / 2 ? j - rows : j);
    const Real kYSqr = kY * kY;

    for (int64_t i = 0; i < cols; ++i)
    {
      const Real kX = deltaKX * (i > cols / 2 ? i - cols : i);
      const Real kXSqr = kX * kX;

      // normalized spatial frequency squared
      const Real radiusSqr = kXSqr / kXNyquistSqr + kYSqr / kYNyquistSqr;

      // band limit the signal to 2/3 of k_Nyquist
      if (9 * radiusSqr > 4)
        result.pixel(i, j) = 0;
      else
        // correct the FFT-scaling
        result.pixel(i, j) = 1 / numPixels;
    }
  }

  return result;
}

RBuffer2D CBandlimitSmooth::mask(const Buffer2DInfo& info)
{
  // We are using a fermi function with mu = 2/3 - alpha and k_B T = 0.015
  constexpr Real alpha = R(0.05);
  constexpr Real mu = R(2.0) / R(3.0) - alpha;
  // Inverse of k_B T
  constexpr Real factor = R(1.0) / R(0.015);

  const int64_t cols = info.cols();
  const int64_t rows = info.rows();
  const Real deltaKX = info.deltaKX();
  const Real deltaKY = info.deltaKY();
  const Real kXNyquistSqr = powerOf<2>(info.kXNyquist());
  const Real kYNyquistSqr = powerOf<2>(info.kYNyquist());
  const Real numPixels = info.numPixels<Real>();

  RBuffer2D result(cols, rows);

  #pragma omp parallel for
  for (int64_t j = 0; j < rows; ++j)
  {
    // coordinates in reciprocal space
    const Real kY = deltaKY * (j > rows / 2 ? j - rows : j);
    const Real kYSqr = kY * kY;

    for (int64_t i = 0; i < cols; ++i)
    {
      const Real kX = deltaKX * (i > cols / 2 ? i - cols : i);
      const Real kXSqr = kX * kX;

      // normalized spatial frequency
      const Real r = std::sqrt(kXSqr / kXNyquistSqr + kYSqr / kYNyquistSqr);

      // correct the FFT-scaling on the fly
      result.pixel(i, j) =
        R(1.0) / ((1 + std::exp((r - mu) * factor)) * numPixels);
    }
  }

  return result;
}

CBandlimitFft::CBandlimitFft(const Buffer2DInfo& info,
                             const CBuffer2D& bufferIn, CBuffer2D& bufferOut,
                             RBuffer2D&& mask) :
  CBandlimit(info, bufferIn, bufferOut),
  mMask(std::move(mask)),
  mBufferFft(bufferIn.createCompatibleBuffer()),
  mFft(Fft::create(bufferIn, mBufferFft, Direction::forward, false)),
  mFftInv(Fft::create(mBufferFft, bufferOut, Direction::backward, false))
{}

void CBandlimitFft::apply()
{
  mFft->transform();
  mBufferFft *= mMask;
  mFftInv->transform();
}


RBandlimit::RBandlimit(Bandlimit bandlimit, const Buffer2DInfo& info) :
  mTempBuffer(info.cols(), info.rows()),
  mCBandlimit(CBandlimit::create(bandlimit, info, mTempBuffer, mTempBuffer))
{}

RBandlimit::RBandlimit(const RBandlimit& other) :
  mTempBuffer(other.mTempBuffer.createCompatibleBuffer()),
  mCBandlimit(CBandlimit::create(other.mBandlimit, other.mCBandlimit->info(),
                                 mTempBuffer, mTempBuffer))
{}

void RBandlimit::apply(RBuffer2D& realBuffer)
{
  const int64_t cols = mTempBuffer.cols();
  const int64_t rows = mTempBuffer.rows();

  assert((realBuffer.cols() == cols) && (realBuffer.rows() == rows));

  // copy the passed buffer into the internal buffer and extend the real
  // values to complex values
  mTempBuffer.assign(realBuffer);

  mCBandlimit->apply();

  // copy the internal buffer back to the real buffer
  // todo: one for-loop
  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      realBuffer.pixel(x, y) = mTempBuffer.pixel(x, y).real();
}

} // namespace Aurora
