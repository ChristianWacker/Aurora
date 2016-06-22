//--- Aurora/UnitTest/PhysicsTest.cpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Physics.hpp"

#include "AuroraLib/Energy.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>

namespace Aurora
{

TEST(Physics, Constants)
{
  EXPECT_NEAR(R(0.52917721092e-10), SI::a0  , R(1e-17));

  EXPECT_NEAR(R(0.510998910e6) , Phys::mE  , R(0.1));
  EXPECT_NEAR(R(1.973269717e3) , Phys::hBar, R(1e-4));
  EXPECT_NEAR(R(5.526349599e-3), Phys::eps0, R(1e-10));
  EXPECT_NEAR(R(0.52917721092) , Phys::a0  , R(1e-7));
}

TEST(Physics, Energy)
{
  EXPECT_DOUBLE_EQ(1.0, Energy(0.0).lorentzFactor());
  EXPECT_DOUBLE_EQ(2.0, Energy(Phys::mE).lorentzFactor());

  EXPECT_DOUBLE_EQ(0.0, lorentzFactorToSpeed(1.0));
}

} // namespace Aurora
