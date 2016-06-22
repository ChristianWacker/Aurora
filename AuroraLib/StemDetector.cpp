//--- Aurora/AuroraLib/StemDetector.cpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/StemDetector.hpp"

#include "AuroraLib/Fft.hpp"

namespace Aurora
{

void StemDetector::connect(const CBuffer2D& inputBuffer)
{
  mInputBufferFft = inputBuffer.createCompatibleBuffer();
  mFft = Fft::create(inputBuffer, mInputBufferFft, Direction::forward, false);
}

AdfDetector::AdfDetector(Buffer2DInfo info, Real alpha, Real beta) :
  StemDetector(info), mAlpha(alpha), mBeta(beta),
  mDetectorMask(detectorMask(1 / std::sqrt(mInfo.numPixels<Real>())))
{ }

CBuffer2D AdfDetector::detectorMask(Real scaling) const
{
  const int64_t cols = mInfo.cols();
  const int64_t rows = mInfo.rows();
  const Real deltaThetaX = mInfo.deltaKX() / mInfo.wavenumber();
  const Real deltaThetaY = mInfo.deltaKY() / mInfo.wavenumber();

  CBuffer2D detectorMask(cols, rows);

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
  {
    const Real thetaY = deltaThetaY * (y > rows / 2 ? y - rows : y);

    for (int64_t x = 0; x < cols; ++x)
    {
      const Real thetaX = deltaThetaX * (x > cols / 2 ? x - cols : x);
      const Real thetaSqr = powerOf<2>(thetaX) + powerOf<2>(thetaY);

      // only pixels inside the ring will be accepted
      if ((thetaSqr >= mAlpha) && (thetaSqr < mBeta))
        detectorMask.pixel(x, y) = scaling;
      else
        detectorMask.pixel(x, y) = 0;
    }
  }

  return detectorMask;
}

Real AdfDetector::signal()
{
  // go to Fourier space
  mFft->transform();
  // apply the mask
  mInputBufferFft.assign(mDetectorMask * mInputBufferFft);
  // calculate the probabilty
  return mInputBufferFft.absSqrReduce();
}

} // namespace Aurora
