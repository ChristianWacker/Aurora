//--- Aurora/AuroraLib/Quadrature.cpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Quadrature.hpp"

#include <iostream>
#include <vector>

namespace Aurora
{

HermitePolynomial::HermitePolynomial(size_t order) :
  mOrder(order), mRoots(order), mWeights(order)
{
  // To calculate the roots of the Hermite polynomial, their values are
  // guessed and subsequently improved with Newton's algorithm. To prevent an
  // overflow during the calculation, we use a rescaled definition of the
  // polynomials as described in Numerical Recipes.

  // the maximal number of Newton iterations
  const int maxNewtonIterations = 10;

  const Real epsilon = R(1e-12);

  // The first five negative roots of the Airy function multiplied by -2^(-1/3)
  const Real airyRoots[5] = {R(1.85575708149), R(3.24460762400),
                             R(4.38167123929), R(5.38661378079),
                             R(6.30526300659)};

  for (size_t i = 0; i < (mOrder + 1) / 2; ++i)
  {
    Real root;

    // calculate a first guess
    if (i < 3)
    {
      // use 18.16.17 [DLMF]
      const Real x = R(2.0) * mOrder + R(1.0);
      root = std::pow(x, R(-1.0) / R(6.0)) * airyRoots[i] - std::sqrt(x);
    }
    else
    {
      // the distance between roots does not change very fast
      root = R(2.0) * mRoots[i - 1] - mRoots[i - 2];
    }

    Real h1;
    Real h2;

    // Use Newton iteration for refinement
    int j = 0;
    for (; j < maxNewtonIterations; ++j)
    {
      // calculate a rescaled Hermite polynomial
      h1 = std::pow(Math::pi, -R(0.25));
      h2 = 0;

      for (size_t k = 0; k < mOrder; ++k)
      {
        Real h3 = h2;
        h2 = h1;
        h1 = root * std::sqrt(R(2.0) / (k + 1)) * h2 -
             std::sqrt(static_cast<Real>(k) / (k + 1)) * h3;
      }

      // h1 now contains the value of the Hermite polynomial evaluated at the
      // root
      if (std::fabs(h1) < epsilon)
        break;

      // Newton step. The derivative is given by
      // H_n' = \sqrt{2 n} H_{n-1}
      root = root - h1 / (std::sqrt(R(2.0) * mOrder) * h2);
    }

    assert(j <= maxNewtonIterations &&
           "Newton iterations did not reach desired accuracy");

    mRoots[i] = root;
    mRoots[mOrder - i - 1] = -root;

    // calculate weight. In our convention it is given by
    // w_i = 2 / H_n'(x_i)^2
    mWeights[i] = R(1.0) / (mOrder * h2 * h2);
    mWeights[mOrder - i - 1] = mWeights[i];
  }
}

void HermitePolynomial::dump(std::ostream& outStream)
{
  outStream << "Hermite-Polynomial Order: " << mOrder << '\n';
  outStream << "Index Root Weight\n";
  for (size_t i = 0; i < mOrder; ++i)
    outStream << i << ' ' << mRoots[i] << " " << mWeights[i] << '\n';
}

/// Internal class, that implements the midpoint rule of integration and a
/// refinement method.
class MidpointRule
{
public:
  MidpointRule(const std::function<Real (Real)>& integrand, Real a, Real b) :
    mIntegrand(integrand), mA(a), mB(b), mOrder(1)
  {
    mResults.push_back((b - a) * integrand(R(0.5) * (a + b)));
  }

  Real refine()
  {
    ++mOrder;
    // 2 * n is the number of sampling that we will add
    int n = 1;

    for (size_t i = 0; i < mOrder - 1; ++i)
      n *= 3;

    // The distance of two sampling points is
    Real delta = (mB - mA) / R(3.0) / n;
    Real twoDelta = R(2.0) * delta;

    // the position of the first new sampling point
    Real x = mA + R(0.5) * delta;
    Real sum = R(0.0);

    // evaluate the integrand at the new sampling points
    for (int i = 0; i < n; ++i)
    {
      sum += mIntegrand(x);
      x += twoDelta;
      sum += mIntegrand(x);
      x += delta;
    }

    Real result = sum * delta + mResults.back() / R(3.0);
    mResults.push_back(result);
    return result;
  }

  Real result(size_t order)
  {
    assert(order <= mOrder);
    return mResults[order - 1];
  }

  size_t order() const
  {
    return mOrder;
  }

private:
  const std::function<Real (Real)> mIntegrand;
  const Real mA;
  const Real mB;
  size_t mOrder;
  std::vector<Real> mResults;
};

Real integrateOpen(const std::function<Real (Real)>& integrand,
                   Real a, Real b, int maxIterations, Real epsilon)
{
  // This function calls MidpointRule::refine() to calculate an adaptive result.
  MidpointRule midpointRule(integrand, a, b);

  // begin with four refinement steps
  Real oldT = 0.0;
  Real oldS = 0.0;

  // If the value of the integration does not change too much with additional
  // refinement steps then return the result.
  for (int i = 0; i < maxIterations; ++i)
  {
    const Real t = midpointRule.refine();
    const Real s = (R(9.0) * t - oldT) / R(8.0);

    if (i > 4)
    {
      if (std::fabs(s - oldS) < epsilon * std::fabs(oldS) ||
          ((R(0.0) == s) && (R(0.0) == oldS)))
        return s;
    }

    oldS = s;
    oldT = t;
  }

  std::cout << "Integral does not converge." << std::endl;
  return oldS;
 // AURORA_THROW(EMath, "Integral does not converge.");
}

class DERule
{
public:
  DERule(const std::function<Real (Real)>& integrand, Real a, Real b,
         Real hMax) :
    mIntegrand(integrand), mA(a), mB(b), mOrder(1), mHMax(hMax)
  {
    Real factor = R(0.25);
    mResults.push_back(hMax * R(2.0) * (b - a) * factor *
                       integrand(R(0.5) * (a + b)));
  }

  Real refine()
  {
    ++mOrder;

    int it = 1;
    for (size_t i = 1; i < mOrder - 1; ++i)
      it *= 2;

    const Real twoH = mHMax / it;
    Real t = R(0.5) * twoH;

    Real sum = 0.0;
    for (int i = 1; i < it; ++i)
    {
      Real q(std::exp(-R(2.0) * std::sinh(t)));
      Real delta((mB - mA) * q / (R(1.0) + q));
      Real fact(q / powerOf<2>(R(1.0) + q) * std::cosh(t));
      sum += fact * (mIntegrand(mA + delta) + mIntegrand(mB - delta));
      t += twoH;
    }

    Real result = R(0.5) * mResults.back() + (mB - mA) * twoH * sum;
    mResults.push_back(result);
    return result;
  }

  Real result(size_t order)
  {
    assert(order <= mOrder);
    return mResults[order - 1];
  }

  size_t order() const
  {
    return mOrder;
  }

private:
  const std::function<Real (Real)> mIntegrand;
  const Real mA;
  const Real mB;
  size_t mOrder;
  const Real mHMax;
  std::vector<Real> mResults;
};

Real integrateDE(const std::function<Real (Real)>& integrand, Real a, Real b,
                 size_t maxIterations, Real epsilon, Real hMax)
{
  // This function calls DERule::refine() to calculate an adaptive result.
  DERule rule(integrand, a, b, hMax);

  Real result = rule.refine();
  Real oldResult;

  // If the value of the integration does not change too much with additional
  // refinement steps then return the result.
  for (size_t i = 0; i < maxIterations; ++i)
  {
    oldResult = result;
    result = rule.refine();

    if ((i > 2) && std::fabs(result - oldResult) <= epsilon * std::fabs(result))
      return result;
  }

  std::cout << "Integral does not converge." << std::endl;
  return result;
}

} // namespace Aurora
