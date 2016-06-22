//--- Aurora/UnitTest/FftTest.cpp ----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Fft.hpp"

#if (AURORA_USE_FFTW == 1)
# include "AuroraLib/FftFftw.hpp"
#endif

#include "AuroraLib/FftKiss.hpp"
#include "AuroraLib/Formatter.hpp"
#include "UnitTest/RelativeError.hpp"

#include <gtest/gtest.h>

namespace Aurora
{

namespace
{
#if (1 == AURORA_SINGLE_PRECISION)
  const Real epsilon = R(1e-5);
#else
  const Real epsilon = R(1e-14);
#endif
}

template<class T>
void checkForwardBackward(const CBuffer2D& testData,
                          CBuffer2D& inBuffer, CBuffer2D& outBuffer,
                          Real epsilon)
{
  // copy the test data into the input array
  inBuffer.assign(testData);

  // create the FFTs
  auto fftForward =
    std::make_unique<T>(inBuffer, outBuffer, Direction::forward);
  auto fftBackward =
    std::make_unique<T>(outBuffer, inBuffer, Direction::backward);

  // do the transformations
  fftForward->transform();
  fftBackward->transform();

  // After a forward and a backward transformation all elements should be
  // simply scaled by numElements.
  const Real numPixelsReal = static_cast<Real>(testData.numPixels());
  for (int64_t y = 0; y < testData.rows(); ++y)
  {
    for (int64_t x = 0; x < testData.cols(); ++x)
    {
      EXPECT_RELATIVE(testData.pixel(x, y).real() * numPixelsReal,
                      inBuffer.pixel(x, y).real(),
                      epsilon)
        << "Element " << x << ", " << y << " does not match\n";
      EXPECT_RELATIVE(testData.pixel(x, y).imag() * numPixelsReal,
                      inBuffer.pixel(x, y).imag(),
                      epsilon)
        << "Element " << x << ", " << y << " does not match\n";
    }
  }
}

void generateTestData2D(CBuffer2D& testData)
{
  const int64_t cols = testData.cols();
  const int64_t rows = testData.rows();

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
      testData.pixel(x, y) =
        Complex(static_cast<Real>(x) / static_cast<Real>(cols) + R(0.5),
                 y < rows / 2 ? R(1.0) : R(2.0));
    }
  }
}

template<class T>
void additionalFftTest()
{
  // Additional checks of the FFT implementations
  std::vector<uint32_t> resolutions = {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                       10, 11, 17, 20, 31, 60};

  try
  {
    // test two-dimensional FFT
    for (uint32_t cols : resolutions)
    {
      for (uint32_t rows : resolutions)
      {
        CBuffer2D testData(cols, rows);
        generateTestData2D(testData);

        auto inBuffer = testData.createCompatibleBuffer();
        auto outBuffer = testData.createCompatibleBuffer();

        checkForwardBackward<T>(testData, inBuffer, outBuffer, epsilon);
      }
    }
  }
  catch (Exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
    ADD_FAILURE();
  }
}

TEST(Fft, Additional)
{
  std::cout << "  Check Kiss...\n";
  additionalFftTest<FftKiss>();

  #if (1 == AURORA_USE_FFTW)
  std::cout << "  Check FFTW-Estimate...\n";
  additionalFftTest<FftFftw<FFTW_ESTIMATE>>();
  std::cout << "  Check FFTW-Measure...\n";
  additionalFftTest<FftFftw<FFTW_MEASURE>>();
  //std::cout << "  Check FFTW-Patient...\n";
  //additionalFftTest<FftFftw<FFTW_PATIENT>>();
  #endif
}

} // namespace Aurora
