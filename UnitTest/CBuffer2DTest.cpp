//--- Aurora/UnitTest/CBuffer2DTest.cpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CBuffer2D.hpp"

#include "AuroraLib/Math.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <random>

namespace Aurora
{

template<int64_t margin>
void generateTestData(CMarginBuffer2D<margin>& buffer)
{
  const int64_t cols = buffer.cols();
  const int64_t rows = buffer.rows();

  // generate test data
  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
      buffer.pixel(x, y) =
        Complex(std::sin(static_cast<Real>(y)), static_cast<Real>(x));
    }
  }
}

template<int64_t margin>
void testCMarginBuffer()
{
  const int64_t cols     = 23;
  const int64_t rows        = 54;
  const int64_t numPixels   = cols * rows;
  const int64_t numElements = (cols + 2 * margin) * (rows + 2 * margin);

  CMarginBuffer2D<margin> buffer1(cols, rows);

  EXPECT_EQ(cols, buffer1.cols());
  EXPECT_EQ(rows, buffer1.rows());
  EXPECT_EQ(numPixels, buffer1.numPixels());
  EXPECT_EQ(numElements, buffer1.numElements());

  // generate test data
  generateTestData(buffer1);

  // check the data in the buffer
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
      EXPECT_EQ(Complex(std::sin(static_cast<Real>(y)), static_cast<Real>(x)),
                buffer1.pixel(x, y));
    }
  }

  // create a compatible buffer
  auto buffer2 = buffer1.createCompatibleBuffer();

  EXPECT_EQ(cols, buffer1.cols());
  EXPECT_EQ(rows, buffer2.rows());
  EXPECT_EQ(numPixels, buffer2.numPixels());
  EXPECT_EQ(numElements, buffer2.numElements());
  EXPECT_TRUE(buffer1.compatible(buffer2));
  EXPECT_TRUE(buffer2.compatible(buffer1));

  // copy the data in the new buffer
  buffer2.assign(buffer1);

  // clear the first buffer
  const Complex clearValue(42.0, -8.0);
  buffer1.setValue(clearValue);

  // check the data in the first buffer
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(clearValue, buffer1.pixel(x, y));

  // check the data in the second buffer
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
       EXPECT_EQ(Complex(std::sin(static_cast<Real>(y)), static_cast<Real>(x)),
                 buffer2.pixel(x, y));
    }
  }

  // test addAssign
  buffer1.setValue(Complex(1.0, -1.0));
  generateTestData(buffer2);
  buffer1 += buffer2;
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
      EXPECT_EQ(Complex(std::sin(static_cast<Real>(y)) + R(1.0),
                        static_cast<Real>(x) - R(1.0)),
                        buffer1.pixel(x, y));
    }
  }

  // test multiplyAssign(CMarginBuffer)
  buffer1.setValue(Complex(1.0, -1.0));
  generateTestData(buffer2);
  buffer1 *= buffer2;
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
      EXPECT_EQ(Complex(std::sin(static_cast<Real>(y)) + static_cast<Real>(x),
                        -std::sin(static_cast<Real>(y)) + static_cast<Real>(x)),
                        buffer1.pixel(x, y));
     }
  }

  // test multiplyAssign(Complex)
  generateTestData(buffer1);
  buffer1 *= Complex(1.0, -1.0);
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
      EXPECT_EQ(Complex(std::sin(static_cast<Real>(y)) + static_cast<Real>(x),
                        -std::sin(static_cast<Real>(y)) + static_cast<Real>(x)),
                        buffer1.pixel(x, y));
     }
  }

  // test multiplyAssign(Real)
  generateTestData(buffer1);
  buffer1 *= 5.0;
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
      EXPECT_EQ(Complex(std::sin(static_cast<Real>(y)) * R(5.0),
                        static_cast<Real>(x) * R(5.0)),
                        buffer1.pixel(x, y));
     }
  }

  // create a third buffer
  const Complex initialValue(-49.0, 7.0);
  CMarginBuffer2D<margin> buffer3(cols + 1, rows, initialValue);
  EXPECT_FALSE(buffer3.compatible(buffer1));
  EXPECT_FALSE(buffer1.compatible(buffer3));

  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(initialValue, buffer3.pixel(x, y));

  std::mt19937 generator;
  std::uniform_real_distribution<Real> distribution(0.0, Math::twoPi);
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols + 1; ++x)
      buffer3.pixel(x, y) = fromArg(distribution(generator));

  EXPECT_EQ(buffer3.numPixels(), buffer3.absSqrReduce());

  // todo: test multiplyAssign
  // todo: test info
}

TEST(CBuffer2D, General)
{
  testCMarginBuffer<0>();
  testCMarginBuffer<1>();
}

TEST(CBuffer2D, Arithmetic)
{
  const int64_t cols = 16;
  const int64_t rows = 8;

  CBuffer2D buffer1(cols, rows, Complex(3.0,  4.0));
  CBuffer2D buffer2(cols, rows, Complex(3.0, -4.0));

  CBuffer2D bufferSum1(cols, rows);
  bufferSum1.assign(buffer1 + buffer2);
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(Complex(6.0, 0.0), bufferSum1.pixel(x, y));

  CBuffer2D bufferSum2(cols, rows);
  bufferSum2.assign(buffer1);
  bufferSum2 += buffer2;
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(Complex(6.0, 0.0), bufferSum2.pixel(x, y));

  CBuffer2D bufferProduct1(cols, rows);
  bufferProduct1.assign(buffer1 * buffer2);
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(Complex(25.0, 0.0), bufferProduct1.pixel(x, y));

  CBuffer2D bufferProduct2(cols, rows);
  bufferProduct2.assign(buffer1);
  bufferProduct2 *= buffer2;
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(Complex(25.0, 0.0), bufferProduct2.pixel(x, y));

  CBuffer2D bufferScale1(cols, rows);
  bufferScale1.assign(2.0 * buffer1);
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(Complex(6.0, 8.0), bufferScale1.pixel(x, y));

  CBuffer2D bufferScale2(cols, rows);
  bufferScale2.assign(buffer1);
  bufferScale2 *= 0.5;
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(Complex(1.5, 2.0), bufferScale2.pixel(x, y));
}

TEST(CBuffer2D, Margin)
{
  const int64_t cols = 10;
  const int64_t rows = 15;

  CMarginBuffer2D<3> buffer(cols, rows, Complex(0.0));

  // generate test data
  for (int64_t y = 0; y < rows; ++y)
  {
    for (int64_t x = 0; x < cols; ++x)
    {
      buffer.pixel(x, y) = Complex(static_cast<Real>(x + y),
                                   static_cast<Real>(y * cols + x));
    }
  }

  buffer.replicateMargin();

  for (int64_t y = -3; y < rows + 3; ++y)
  {
    EXPECT_EQ(buffer.pixel(cols - 2, y), buffer.pixel(-2, y));
    EXPECT_EQ(buffer.pixel(2, y), buffer.pixel(cols + 2, y));
  }

  for (int64_t x = -3; x < cols + 3; ++x)
  {
    EXPECT_EQ(buffer.pixel(x, rows - 2), buffer.pixel(x, -2));
    EXPECT_EQ(buffer.pixel(x, 2), buffer.pixel(x, rows + 2));
  }
}

TEST(CBuffer2D, FftShiftOdd)
{
  auto original = std::make_unique<CBuffer2D>(3, 3);
  auto originalData = original->data();

  // generate test data
  for (int64_t i = 0; i < 9; ++i)
    originalData[i] = Complex(static_cast<Real>(i));

  auto shifted = original->fftShift();
  auto shiftedData = shifted.data();

  EXPECT_EQ(Complex(8.0), shiftedData[0]);
  EXPECT_EQ(Complex(6.0), shiftedData[1]);
  EXPECT_EQ(Complex(7.0), shiftedData[2]);
  EXPECT_EQ(Complex(2.0), shiftedData[3]);
  EXPECT_EQ(Complex(0.0), shiftedData[4]);
  EXPECT_EQ(Complex(1.0), shiftedData[5]);
  EXPECT_EQ(Complex(5.0), shiftedData[6]);
  EXPECT_EQ(Complex(3.0), shiftedData[7]);
  EXPECT_EQ(Complex(4.0), shiftedData[8]);
}

TEST(CBuffer2D, FftShiftEven)
{
  auto original = std::make_shared<CBuffer2D>(4, 4);
  auto originalData = original->data();

  // generate test data
  for (int64_t i = 0; i < 16; ++i)
    originalData[i] = Complex(static_cast<Real>(i));

  auto shifted = original->fftShift();
  auto shiftedData = shifted.data();

  EXPECT_EQ(Complex(10.0), shiftedData[ 0]);
  EXPECT_EQ(Complex(11.0), shiftedData[ 1]);
  EXPECT_EQ(Complex( 8.0), shiftedData[ 2]);
  EXPECT_EQ(Complex( 9.0), shiftedData[ 3]);
  EXPECT_EQ(Complex(14.0), shiftedData[ 4]);
  EXPECT_EQ(Complex(15.0), shiftedData[ 5]);
  EXPECT_EQ(Complex(12.0), shiftedData[ 6]);
  EXPECT_EQ(Complex(13.0), shiftedData[ 7]);
  EXPECT_EQ(Complex( 2.0), shiftedData[ 8]);
  EXPECT_EQ(Complex( 3.0), shiftedData[ 9]);
  EXPECT_EQ(Complex( 0.0), shiftedData[10]);
  EXPECT_EQ(Complex( 1.0), shiftedData[11]);
  EXPECT_EQ(Complex( 6.0), shiftedData[12]);
  EXPECT_EQ(Complex( 7.0), shiftedData[13]);
  EXPECT_EQ(Complex( 4.0), shiftedData[14]);
  EXPECT_EQ(Complex( 5.0), shiftedData[15]);
}


} // namespace Aurora
