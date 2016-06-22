//--- Aurora/UnitTest/MathTest.cpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Math.hpp"
#include "UnitTest/RelativeError.hpp"

#include <gtest/gtest.h>

namespace Aurora
{

TEST(math, IsPowerOfTwo)
{
  // test isPowerOfTwo()
  EXPECT_FALSE(isPowerOfTwo(0));
  EXPECT_TRUE (isPowerOfTwo(1));
  EXPECT_TRUE (isPowerOfTwo(2));
  EXPECT_FALSE(isPowerOfTwo(3));
  EXPECT_TRUE (isPowerOfTwo(4));
  EXPECT_FALSE(isPowerOfTwo(65535));
  EXPECT_TRUE (isPowerOfTwo(65536));
  EXPECT_FALSE(isPowerOfTwo(65537));
  EXPECT_FALSE(isPowerOfTwo(65538));
}

TEST(Math, PowerOf)
{
  // test power_of template
  EXPECT_DOUBLE_EQ(std::pow(0.3, 1), powerOf<1>(0.3));
  EXPECT_DOUBLE_EQ(std::pow(0.3, 2), powerOf<2>(0.3));
  EXPECT_DOUBLE_EQ(std::pow(0.3, 4), powerOf<4>(0.3));
  EXPECT_DOUBLE_EQ(std::pow(0.3, 9), powerOf<9>(0.3));

  EXPECT_FLOAT_EQ(std::pow(0.3f, 1.0f), powerOf<1>(0.3f));
  EXPECT_FLOAT_EQ(std::pow(0.3f, 2.0f), powerOf<2>(0.3f));
  EXPECT_FLOAT_EQ(std::pow(0.3f, 4.0f), powerOf<4>(0.3f));
  EXPECT_FLOAT_EQ(std::pow(0.3f, 9.0f), powerOf<9>(0.3f));
}

TEST(Math, Factorial)
{
  EXPECT_EQ(1, factorial(0));
  EXPECT_EQ(1, factorial(1));
  EXPECT_EQ(2, factorial(2));
  EXPECT_EQ(6, factorial(3));
  EXPECT_EQ(720, factorial(6));
  EXPECT_EQ(1307674368000, factorial(15));

}

TEST(Math, Binomial)
{
  const Real epsilon = R(1e-7);
  EXPECT_RELATIVE(R(  0.0       ), binomial(R( 1.0), -1), epsilon);
  EXPECT_RELATIVE(R(  1.0       ), binomial(R(10.0),  0), epsilon);
  EXPECT_RELATIVE(R( 10.0       ), binomial(R(10.0),  1), epsilon);
  EXPECT_RELATIVE(R( 45.0       ), binomial(R(10.0),  2), epsilon);
  EXPECT_RELATIVE(R(252.0       ), binomial(R(10.0),  5), epsilon);
  EXPECT_RELATIVE(R(  0.02734375), binomial(R( 0.5),  5), epsilon);
}

TEST(Math, BesselI0)
{
  const Real epsilon = R(1e-6);

  // test values from Mathematica 9
  EXPECT_RELATIVE(R(2.81571662847e3), besselI0(R(-10.0)), epsilon);
  EXPECT_RELATIVE(R(1.06348337074  ), besselI0(R( -0.5)), epsilon);
  EXPECT_RELATIVE(R(1.00250156293  ), besselI0(R( -0.1)), epsilon);
  EXPECT_RELATIVE(R(1.0            ), besselI0(R(  0.0)), epsilon);
  EXPECT_RELATIVE(R(1.00250156293  ), besselI0(R(  0.1)), epsilon);
  EXPECT_RELATIVE(R(1.06348337074  ), besselI0(R(  0.5)), epsilon);
  EXPECT_RELATIVE(R(2.81571662847e3), besselI0(R( 10.0)), epsilon);
}

TEST(Math, BesselI1)
{
  const Real epsilon = R(1e-6);

  // test values from Mathematica 9
  EXPECT_RELATIVE(R(-2.67098830370e3 ), besselI1(R(-10.0)), epsilon);
  EXPECT_RELATIVE(R(-2.57894305391e-1), besselI1(R( -0.5)), epsilon);
  EXPECT_RELATIVE(R(-5.00625260471e-2), besselI1(R( -0.1)), epsilon);
  EXPECT_RELATIVE(R( 0.0             ), besselI1(R(  0.0)), epsilon);
  EXPECT_RELATIVE(R( 5.00625260471e-2), besselI1(R(  0.1)), epsilon);
  EXPECT_RELATIVE(R( 2.57894305391e-1), besselI1(R(  0.5)), epsilon);
  EXPECT_RELATIVE(R( 2.67098830370e3 ), besselI1(R( 10.0)), epsilon);
}

TEST(Math, BesselK0)
{
  const Real epsilon = R(1e-6);

  // test values from Mathematica 9
  EXPECT_TRUE(std::isinf(besselK0(R(0.0))));
  EXPECT_RELATIVE(R(2.42706902470   ), besselK0(R( 0.1)), epsilon);
  EXPECT_RELATIVE(R(0.924419071228  ), besselK0(R( 0.5)), epsilon);
  EXPECT_RELATIVE(R(1.77800623162e-5), besselK0(R(10.0)), epsilon);
}

TEST(Math, BesselK1)
{
  const Real epsilon = R(1e-6);

  // test values from Mathematica 9
  EXPECT_TRUE(std::isinf(besselK1(R(0.0))));
  EXPECT_RELATIVE(R(9.85384478087   ), besselK1(R( 0.1)), epsilon);
  EXPECT_RELATIVE(R(1.65644112000   ), besselK1(R( 0.5)), epsilon);
  EXPECT_RELATIVE(R(1.86487734538e-5), besselK1(R(10.0)), epsilon);
}

TEST(Vector3R, ArithmeticOperators)
{
  Vector3R vec1(R(0.5 ), R(1.0), R(2.0));
  Vector3R vec2(R(0.25), R(1.0), R(4.0));

  EXPECT_EQ(vec1, vec2 / vec1);
}

TEST(Vector3R, MinimumAndMaximum)
{
  Vector3R vec1(R(0.5), R( 1.0), R(2.0));
  Vector3R vec2(R(2.0), R(-1.0), R(0.5));

  EXPECT_EQ(Vector3R(R(2.0), R( 1.0), R(2.0)), max(vec1, vec2));
  EXPECT_EQ(Vector3R(R(0.5), R(-1.0), R(0.5)), min(vec1, vec2));
}

} // namespace Aurora
