//--- Aurora/UnitTest/StatisticsTest.cpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Statistics.hpp"

#include <gtest/gtest.h>

namespace Aurora
{

TEST(Statistics, General)
{
  Statistics<Real> stats;
  stats.add(9.0);
  stats.add(10.0);
  stats.add(14.0);

  EXPECT_EQ(11.0, stats.mean());
  EXPECT_EQ(10.0, stats.median());
  EXPECT_EQ(9.0, stats.min());
  EXPECT_EQ(14.0, stats.max());
  EXPECT_EQ(7.0, stats.variance());
}

} // namespace Aurora
