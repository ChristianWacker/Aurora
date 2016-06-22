//--- Aurora/UnitTest/UtilsTest.cpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Utils.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>

namespace Aurora
{

TEST(Utils, clamp)
{
  // test clamp()
  EXPECT_EQ( 0, clamp<int>(-1, 0 , 10));
  EXPECT_EQ( 5, clamp<int>( 5, 0 , 10));
  EXPECT_EQ(10, clamp<int>(11, 0 , 10));
}

TEST(Utils, saturate)
{
  EXPECT_EQ(0, saturate(std::numeric_limits<int>::min()));
  EXPECT_EQ(0, saturate(-1));
  EXPECT_EQ(0, saturate(0));
  EXPECT_EQ(1, saturate(1));
  EXPECT_EQ(128, saturate(128));
  EXPECT_EQ(255, saturate(255));
  EXPECT_EQ(255, saturate(256));
  EXPECT_EQ(255, saturate(std::numeric_limits<int>::max()));
}

TEST(Utils, trim)
{
  std::string reference("reference");
  EXPECT_EQ(reference, trim(" \t reference"));
  EXPECT_EQ(reference, trim(" reference\t"));
}

} // namespace Aurora
