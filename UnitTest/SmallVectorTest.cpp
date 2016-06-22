//--- Aurora/UnitTest/SmallVectorTest.cpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/SmallVector.hpp"

#include <gtest/gtest.h>

namespace Aurora
{

TEST(SmallVector, Construction)
{
  SmallVector<int, 40> a(30);
  EXPECT_EQ(30, a.size());
  SmallVector<int, 40> b(40);
  EXPECT_EQ(40, b.size());
  SmallVector<int, 40> c(50);
  EXPECT_EQ(50, c.size());
}

} // namespace Aurora
