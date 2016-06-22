//--- Aurora/UnitTest/RBuffer2DTest.cpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/RBuffer2D.hpp"
#include "AuroraLib/Math.hpp"

#include <cmath>
#include <gtest/gtest.h>

namespace Aurora
{

void generateTestData(RBuffer2D& buffer)
{
  const int64_t cols = buffer.cols();
  const int64_t rows = buffer.rows();

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      buffer.pixel(x, y) = std::sin(static_cast<Real>(y));
}

TEST(RBuffer2D, General)
{
  const uint32_t cols      = 23;
  const uint32_t rows      = 54;
  const uint32_t numPixels = cols * rows;

  RBuffer2D buffer1(cols, rows, 0.0);

  EXPECT_EQ(cols, buffer1.cols());
  EXPECT_EQ(rows, buffer1.rows());
  EXPECT_EQ(numPixels, buffer1.numPixels());
  EXPECT_EQ(numPixels, buffer1.numElements());

  // generate test data
  generateTestData(buffer1);

  // check the data in the buffer
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(std::sin(static_cast<Real>(y)), buffer1.pixel(x, y));

  // create a compatible buffer
  auto buffer2 = buffer1.createCompatibleBuffer();

  EXPECT_EQ(cols, buffer1.cols());
  EXPECT_EQ(rows, buffer2.rows());
  EXPECT_EQ(numPixels, buffer2.numPixels());
  EXPECT_EQ(numPixels, buffer2.numElements());
  EXPECT_TRUE(buffer1.compatible(buffer2));
  EXPECT_TRUE(buffer2.compatible(buffer1));

  // copy the data in the new buffe
  buffer2.assign(buffer1);

  // clear the first buffer
  Real clearValue = 42.0;
  buffer1.setValue(clearValue);

  // check the data in the first buffer
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(clearValue, buffer1.pixel(x, y));

  // check the data in the second buffer
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(std::sin(static_cast<Real>(y)), buffer2.pixel(x, y));

  // test addAssign
  buffer1.setValue(-1.0);
  generateTestData(buffer2);
  buffer1 += buffer2;
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(std::sin(static_cast<Real>(y)) - R(1.0), buffer1.pixel(x, y));

  // test multiplyAssign(RBuffer2D)
  buffer1.setValue(-1.0);
  generateTestData(buffer2);
  buffer1 *= buffer2;
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(-std::sin(static_cast<Real>(y)), buffer1.pixel(x, y));

  // test multiplyAssign(RealD)
  generateTestData(buffer1);
  buffer1 *= -1.0;
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(-std::sin(static_cast<Real>(y)), buffer1.pixel(x, y));

  // create a third buffer
  const Real initialValue = -49.0;
  RBuffer2D buffer3(cols + 1, rows, initialValue);
  EXPECT_FALSE(buffer3.compatible(buffer1));
  EXPECT_FALSE(buffer1.compatible(buffer3));

  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      EXPECT_EQ(initialValue, buffer3.pixel(x, y));

  // todo: test info
}

} // namespace Aurora
