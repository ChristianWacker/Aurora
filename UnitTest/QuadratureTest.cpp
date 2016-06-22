//--- Aurora/UnitTest/QuadratureTest.cpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Quadrature.hpp"
#include "UnitTest/RelativeError.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>

namespace Aurora
{

TEST(Quadrature, Hermite)
{
#if (1 == AURORA_SINGLE_PRECISION)
  const Real epsilon = R(1e-6);
#else
  const Real epsilon = R(1e-9);
#endif

  EXPECT_RELATIVE(R(  -4.0)        , hermite(3, R( 1.0)), epsilon);
  EXPECT_RELATIVE(R(2141.371189760), hermite(8, R(-1.3)), epsilon);
}

TEST(Quadrature, Roots)
{
#if (1 == AURORA_SINGLE_PRECISION)
  const Real epsilon = R(1e-6);
#else
  const Real epsilon = R(1e-9);
#endif

  // values of the roots calculated with the website
  // http://www.efunda.com/math/num_integration/findgausshermite.cfm

  HermitePolynomial hermiteInfo4(4);
  EXPECT_RELATIVE(R(-1.65068012389 ), hermiteInfo4.root(0), epsilon);
  EXPECT_RELATIVE(R(-0.524647623275), hermiteInfo4.root(1), epsilon);
  EXPECT_RELATIVE(R( 0.524647623275), hermiteInfo4.root(2), epsilon);
  EXPECT_RELATIVE(R( 1.65068012389 ), hermiteInfo4.root(3), epsilon);

  HermitePolynomial hermiteInfo5(5);
  EXPECT_RELATIVE(R(-2.02018287046 ), hermiteInfo5.root(0), epsilon);
  EXPECT_RELATIVE(R(-0.958572464614), hermiteInfo5.root(1), epsilon);
  EXPECT_RELATIVE(R( 0.0           ), hermiteInfo5.root(2), epsilon);
  EXPECT_RELATIVE(R( 0.958572464614), hermiteInfo5.root(3), epsilon);
  EXPECT_RELATIVE(R( 2.02018287046 ), hermiteInfo5.root(4), epsilon);

  HermitePolynomial hermiteInfo7(7);
  EXPECT_RELATIVE(R(-2.65196135684 ), hermiteInfo7.root(0), epsilon);
  EXPECT_RELATIVE(R(-1.67355162877 ), hermiteInfo7.root(1), epsilon);
  EXPECT_RELATIVE(R(-0.816287882859), hermiteInfo7.root(2), epsilon);
  EXPECT_RELATIVE(R( 0.0           ), hermiteInfo7.root(3), epsilon);
  EXPECT_RELATIVE(R( 0.816287882859), hermiteInfo7.root(4), epsilon);
  EXPECT_RELATIVE(R( 1.67355162877 ), hermiteInfo7.root(5), epsilon);
  EXPECT_RELATIVE(R( 2.65196135684 ), hermiteInfo7.root(6), epsilon);

  HermitePolynomial hermiteInfo32(32);
  EXPECT_RELATIVE(R(-7.12581390983 ), hermiteInfo32.root( 0), epsilon);
  EXPECT_RELATIVE(R(-6.40949814928 ), hermiteInfo32.root( 1), epsilon);
  EXPECT_RELATIVE(R(-5.81222594946 ), hermiteInfo32.root( 2), epsilon);
  EXPECT_RELATIVE(R(-5.27555098664 ), hermiteInfo32.root( 3), epsilon);
  EXPECT_RELATIVE(R(-4.77716450334 ), hermiteInfo32.root( 4), epsilon);
  EXPECT_RELATIVE(R(-4.30554795347 ), hermiteInfo32.root( 5), epsilon);
  EXPECT_RELATIVE(R(-3.85375548542 ), hermiteInfo32.root( 6), epsilon);
  EXPECT_RELATIVE(R(-3.41716749282 ), hermiteInfo32.root( 7), epsilon);
  EXPECT_RELATIVE(R(-2.99249082501 ), hermiteInfo32.root( 8), epsilon);
  EXPECT_RELATIVE(R(-2.57724953773 ), hermiteInfo32.root( 9), epsilon);
  EXPECT_RELATIVE(R(-2.16949918361 ), hermiteInfo32.root(10), epsilon);
  EXPECT_RELATIVE(R(-1.76765410946 ), hermiteInfo32.root(11), epsilon);
  EXPECT_RELATIVE(R(-1.37037641095 ), hermiteInfo32.root(12), epsilon);
  EXPECT_RELATIVE(R(-0.97650046359 ), hermiteInfo32.root(13), epsilon);
  EXPECT_RELATIVE(R(-0.584978765436), hermiteInfo32.root(14), epsilon);
  EXPECT_RELATIVE(R(-0.194840741569), hermiteInfo32.root(15), epsilon);
  EXPECT_RELATIVE(R( 0.194840741569), hermiteInfo32.root(16), epsilon);
  EXPECT_RELATIVE(R( 0.584978765436), hermiteInfo32.root(17), epsilon);
  EXPECT_RELATIVE(R( 0.97650046359 ), hermiteInfo32.root(18), epsilon);
  EXPECT_RELATIVE(R( 1.37037641095 ), hermiteInfo32.root(19), epsilon);
  EXPECT_RELATIVE(R( 1.76765410946 ), hermiteInfo32.root(20), epsilon);
  EXPECT_RELATIVE(R( 2.16949918361 ), hermiteInfo32.root(21), epsilon);
  EXPECT_RELATIVE(R( 2.57724953773 ), hermiteInfo32.root(22), epsilon);
  EXPECT_RELATIVE(R( 2.99249082501 ), hermiteInfo32.root(23), epsilon);
  EXPECT_RELATIVE(R( 3.41716749282 ), hermiteInfo32.root(24), epsilon);
  EXPECT_RELATIVE(R( 3.85375548542 ), hermiteInfo32.root(25), epsilon);
  EXPECT_RELATIVE(R( 4.30554795347 ), hermiteInfo32.root(26), epsilon);
  EXPECT_RELATIVE(R( 4.77716450334 ), hermiteInfo32.root(27), epsilon);
  EXPECT_RELATIVE(R( 5.27555098664 ), hermiteInfo32.root(28), epsilon);
  EXPECT_RELATIVE(R( 5.81222594946 ), hermiteInfo32.root(29), epsilon);
  EXPECT_RELATIVE(R( 6.40949814928 ), hermiteInfo32.root(30), epsilon);
  EXPECT_RELATIVE(R( 7.12581390983 ), hermiteInfo32.root(31), epsilon);
}

TEST(Quadrature, Weights)
{
#if (1 == AURORA_SINGLE_PRECISION)
  const Real epsilon = R(1e-5);
#else
  const Real epsilon = R(1e-8);
#endif

  // values of the weights calculated with the website
  // http://www.efunda.com/math/num_integration/findgausshermite.cfm

  HermitePolynomial hermiteInfo3(3);
  EXPECT_RELATIVE(R(0.295408975151), hermiteInfo3.weight(0), epsilon);
  EXPECT_RELATIVE(R(1.1816359006  ), hermiteInfo3.weight(1), epsilon);
  EXPECT_RELATIVE(R(0.295408975151), hermiteInfo3.weight(2), epsilon);

  HermitePolynomial hermiteInfo6(6);
  EXPECT_RELATIVE(R(0.00453000990551), hermiteInfo6.weight(0), epsilon);
  EXPECT_RELATIVE(R(0.157067320323  ), hermiteInfo6.weight(1), epsilon);
  EXPECT_RELATIVE(R(0.724629595224  ), hermiteInfo6.weight(2), epsilon);
  EXPECT_RELATIVE(R(0.724629595224  ), hermiteInfo6.weight(3), epsilon);
  EXPECT_RELATIVE(R(0.157067320323  ), hermiteInfo6.weight(4), epsilon);
  EXPECT_RELATIVE(R(0.00453000990551), hermiteInfo6.weight(5), epsilon);

  HermitePolynomial hermiteInfo31(31);
  EXPECT_RELATIVE(R(4.61896839487e-22), hermiteInfo31.weight( 0), epsilon);
  EXPECT_RELATIVE(R(5.11060900411e-18), hermiteInfo31.weight( 1), epsilon);
  EXPECT_RELATIVE(R(5.89955650839e-15), hermiteInfo31.weight( 2), epsilon);
  EXPECT_RELATIVE(R(1.86037351576e-12), hermiteInfo31.weight( 3), epsilon);
  EXPECT_RELATIVE(R(2.35249200949e-10), hermiteInfo31.weight( 4), epsilon);
  EXPECT_RELATIVE(R(1.4611988324e-8  ), hermiteInfo31.weight( 5), epsilon);
  EXPECT_RELATIVE(R(5.04371256258e-7 ), hermiteInfo31.weight( 6), epsilon);
  EXPECT_RELATIVE(R(1.0498602756e-5  ), hermiteInfo31.weight( 7), epsilon);
  EXPECT_RELATIVE(R(0.000139520903955), hermiteInfo31.weight( 8), epsilon);
  EXPECT_RELATIVE(R(0.00123368330734 ), hermiteInfo31.weight( 9), epsilon);
  EXPECT_RELATIVE(R(0.0074827999141  ), hermiteInfo31.weight(10), epsilon);
  EXPECT_RELATIVE(R(0.0318472307313  ), hermiteInfo31.weight(11), epsilon);
  EXPECT_RELATIVE(R(0.0967179481609  ), hermiteInfo31.weight(12), epsilon);
  EXPECT_RELATIVE(R(0.212132788669   ), hermiteInfo31.weight(13), epsilon);
  EXPECT_RELATIVE(R(0.338772657894   ), hermiteInfo31.weight(14), epsilon);
  EXPECT_RELATIVE(R(0.395778556099   ), hermiteInfo31.weight(15), epsilon);
  EXPECT_RELATIVE(R(0.338772657894   ), hermiteInfo31.weight(16), epsilon);
  EXPECT_RELATIVE(R(0.212132788669   ), hermiteInfo31.weight(17), epsilon);
  EXPECT_RELATIVE(R(0.0967179481609  ), hermiteInfo31.weight(18), epsilon);
  EXPECT_RELATIVE(R(0.0318472307313  ), hermiteInfo31.weight(19), epsilon);
  EXPECT_RELATIVE(R(0.0074827999141  ), hermiteInfo31.weight(20), epsilon);
  EXPECT_RELATIVE(R(0.00123368330734 ), hermiteInfo31.weight(21), epsilon);
  EXPECT_RELATIVE(R(0.000139520903955), hermiteInfo31.weight(22), epsilon);
  EXPECT_RELATIVE(R(1.0498602756e-5  ), hermiteInfo31.weight(23), epsilon);
  EXPECT_RELATIVE(R(5.04371256258e-7 ), hermiteInfo31.weight(24), epsilon);
  EXPECT_RELATIVE(R(1.4611988324e-8  ), hermiteInfo31.weight(25), epsilon);
  EXPECT_RELATIVE(R(2.35249200949e-10), hermiteInfo31.weight(26), epsilon);
  EXPECT_RELATIVE(R(1.86037351576e-12), hermiteInfo31.weight(27), epsilon);
  EXPECT_RELATIVE(R(5.89955650839e-15), hermiteInfo31.weight(28), epsilon);
  EXPECT_RELATIVE(R(5.11060900411e-18), hermiteInfo31.weight(29), epsilon);
  EXPECT_RELATIVE(R(4.61896839487e-22), hermiteInfo31.weight(30), epsilon);
}

TEST(Quadrature, GaussHermite)
{
#if (1 == AURORA_SINGLE_PRECISION)
  const Real epsilon = R(1e-6);
#else
  const Real epsilon = R(1e-7);
#endif

  // test values calculated with mathematica
  HermitePolynomial hermite4(4);
  HermitePolynomial hermite5(5);
  HermitePolynomial hermite30(30);

  auto f1 = [](Real x) { return x * x; };
  EXPECT_RELATIVE(R(33.0153384200), integrateGaussHermite<Real>(f1, hermite4 , R(5.6789), R(0.7654321)), epsilon);
  EXPECT_RELATIVE(R(33.0153384200), integrateGaussHermite<Real>(f1, hermite5 , R(5.6789), R(0.7654321)), epsilon);
  EXPECT_RELATIVE(R(33.0153384200), integrateGaussHermite<Real>(f1, hermite30, R(5.6789), R(0.7654321)), epsilon);

  auto f2 = [](Real x) { return fromArg(x); };
  Complex result2 = integrateGaussHermite<Complex>(f2, hermite4, R(4.0), R(0.0625));
  EXPECT_RELATIVE(R(-0.6335331208), result2.real(), epsilon);
  EXPECT_RELATIVE(R(-0.7335181304), result2.imag(), epsilon);

//  EXPECT_NEAR(result2, gaussHermiteQuadrature<Complex>(hermite30, 4.0, 0.25, func2), 1e-5);
}

TEST(Quadrature, Integrate)
{
#if (1 == AURORA_SINGLE_PRECISION)
  const Real epsilon = R(1e-4);
#else
  const Real epsilon = R(1e-6);
#endif

  auto f1 = [](Real x) { return std::exp(-x) / std::sqrt(x); };
  Real result1 = R(1.772453850);
  EXPECT_RELATIVE( result1, integrateDE(f1, R( 0.0), R(20.0), 20, epsilon), epsilon);
  EXPECT_RELATIVE(-result1, integrateDE(f1, R(20.0), R( 0.0), 20, epsilon), epsilon);

  auto f2 = [](Real x) { return x * x; };
  Real result2 = R(8000.0) / R(3.0);
  EXPECT_RELATIVE(result2, integrateOpen(f2, R(0.0), R(20.0), 20, epsilon), epsilon);
  EXPECT_RELATIVE(result2, integrateDE(f2, R(0.0), R(20.0)  , 20, epsilon), epsilon);
}

} // namespace Aurora
