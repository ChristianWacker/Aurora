//--- Aurora/UnitTest/HistogramTest.cpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Histogram.hpp"

#include <gtest/gtest.h>

namespace Aurora
{

TEST(Histogram, General)
{
  Histogram hist1(-1.0, 1.0, 10);

  EXPECT_EQ(10, hist1.numBins());
  EXPECT_EQ(-1.0, hist1.minValue());
  EXPECT_EQ(1.0, hist1.maxValue());

  EXPECT_EQ(7, hist1.valueToBin(hist1.binToValue(7)));
  EXPECT_EQ(8, hist1.valueToBin(R(0.79)));
  EXPECT_EQ(9, hist1.valueToBin(R(0.81)));
  EXPECT_EQ(9, hist1.valueToBin(R(0.99)));
  EXPECT_EQ(9, hist1.valueToBin(R(1.0)));


  for (size_t i = 0; i < 10; ++i)
    EXPECT_EQ(0, hist1[i]);

  hist1.insert(R(1e-6));
  EXPECT_EQ(1, hist1[5]);
  EXPECT_EQ(1, hist1.numEntries());

  hist1.insert(R( 0.1));
  hist1.insert(R(-0.1));

  EXPECT_EQ(3, hist1.numEntries());
  EXPECT_EQ(5, hist1.mode());
  EXPECT_DOUBLE_EQ(0.1 / 3.0, hist1.mean());
}

} // namespace Aurora
