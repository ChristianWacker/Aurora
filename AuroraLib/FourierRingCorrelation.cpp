//--- Aurora/AuroraLib/FourierRingCorrelation.cpp ------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/FourierRingCorrelation.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Math.hpp"

namespace Aurora
{

FourierRingCorrelation::FourierRingCorrelation(const CBuffer2D& buffer1,
                                               const CBuffer2D& buffer2,
                                               Real deltaKX, Real deltaKY,
                                               Real binSize)
{
  assert(buffer1.compatible(buffer2));

  // calculate the Fourier transformation
  auto buffer1Fft = buffer1.fft(Direction::forward);
  auto buffer2Fft = buffer2.fft(Direction::forward);

  const int64_t cols = buffer1Fft.cols();
  const int64_t rows = buffer1Fft.rows();

  // the maximal radius is half the diagonal of the rectangle
  const Real maxRadius = std::hypot(R(0.5) * cols * deltaKX,
                                    R(0.5) * rows * deltaKY);
  mNumBins = static_cast<size_t>(maxRadius / binSize) + 1;

  std::vector<Real> a(mNumBins, 0);
  std::vector<Real> b(mNumBins, 0);
  std::vector<Complex> c(mNumBins, Complex(0));

  // Because of the Friedel Symmetry of the FFT of a real valued source, we only
  // need to process one half of the pixels
  for (int64_t kY = 0; kY < rows / 2; ++kY)
  {
    const Real kYReal = deltaKY * kY;

    for (int64_t kX = 0; kX < cols; ++kX)
    {
      // this calculation is needed to shift the lower frequencies of the
      // Fourier transform into the center
      const Real kXReal = deltaKX * ((kX > cols / 2) ? (kX - cols) : kX);

      const Real radius = std::hypot(kXReal, kYReal);

      const size_t index = static_cast<size_t>(radius / binSize);

      c[index] += buffer1Fft.pixel(kX, kY) * conj(buffer2Fft.pixel(kX, kY));
      a[index] += absSqr(buffer1Fft.pixel(kX, kY));
      b[index] += absSqr(buffer2Fft.pixel(kX, kY));
    }
  }

  for (size_t i = 0; i < mNumBins; ++i)
    mResult.push_back(c[i].real() / std::sqrt(a[i] * b[i]));
}

} // namespace Aurora
