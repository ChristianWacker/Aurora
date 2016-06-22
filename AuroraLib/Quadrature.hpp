//--- Aurora/AuroraLib/Quadrature.hpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_QUADRATURE_HPP
#define AURORA_AURORA_LIB_QUADRATURE_HPP

#include "AuroraLib/Math.hpp"

#include <cmath>
#include <functional>
#include <vector>

namespace Aurora
{

/// Calculates the physical Hermite polynomial, defined as
///  @f[
///     H_n(x)=(-1)^n \exp(x^2) \frac{\text{d}^n}{\text{d}x^n} \exp(-x^2)
///  @f]
inline Real hermite(size_t n, Real x)
{
  Real h1 = R(1.0);
  Real h2 = R(0.0);

  for (size_t i = 0; i < n; ++i)
  {
    Real h3 = h2;
    h2 = h1;
    h1 = R(2.0) * (x * h2 - i * h3);
  }

  return h1;
}

/// This class provides information about Hermite polynomials needed by the
/// Gauss-Hermite-Integration algorithm.
class AURORA_API HermitePolynomial
{
public:
  /// Constructur.
  /// @param order
  ///   The order of the Hermite polynomial.
  HermitePolynomial(size_t order);

  /// Returns the order of the polynomial.
  size_t order() const
  {
    return mOrder;
  }

  /// Return the i-th root @f$x_i@f$ of the n-th polynomial
  ///  @f[
  ///      -\infty < x_0 < x_1 < \ldots < x_{n-1} < \infty
  ///  @f]
  Real root(size_t i) const
  {
    assert(i < mOrder);
    return mRoots[i];
  }

  /// Returns the weight @f$w_i@f$ associated with the i-th root @f$x_i@f$.
  Real weight(size_t i) const
  {
    assert(i < mOrder);
    return mWeights[i];
  }

  void dump(std::ostream& outStream);

private:
  const size_t mOrder;
  std::vector<Real> mRoots;
  std::vector<Real> mWeights;
};

/// This function evaluates the integral
///  @f[
///    \frac{1}{\sigma\sqrt{2\,\pi}}\intop_{-\infty}^{\infty} \
///      \exp\left(-\frac{(x-\mu)^2}{2\,\sigma^2}\right) f(x) \text{d}x
///  @f]
/// using the Gauss-Hermite-Integration algorithm. Beside the function @p f
/// a Hermite polynomial @p polynomial has to be specified. This polynomial
/// determines the accuracy of the integration.
template<class T>
T integrateGaussHermite(std::function<T (Real)> f,
                        const HermitePolynomial& polynomial, Real mean,
                        Real variance)
{
  T result = 0;

  const size_t iterations = polynomial.order() / 2;
  for (size_t i = 0; i < iterations; ++i)
  {
    const Real abscissa = std::sqrt(2 * variance) * polynomial.root(i);
    result += polynomial.weight(i) *
              (f(abscissa + mean) + f(-abscissa + mean));
  }

  if (1 == polynomial.order() % 2)
  {
    // odd order polynomial
    result += polynomial.weight(polynomial.order() / 2) * f(mean);
  }

  return result / std::sqrt(Math::pi);
}

/// Integrates the functor @p integrand over the open interval (a, b).
/// @param integrand
///  The integrand that should be used for the integration.
/// @param a
///  The left endpoint of the integration interval.
/// @param b
///  The right endpoint of the integration interval.
/// @param maxIterations
///  The maximal number of iterations that the functions should use to achieve
///  the require precision. In each iteration the number of sampling points is
///  tripled.
/// @param epsilon
///  The upper bound on the integration error. If the function is not able to
///  achieved the required precision, an exception is thrown.
AURORA_API Real integrateOpen(const std::function<Real (Real)>& integrand,
                              Real a, Real b, int maxIterations = 15,
                              Real epsilon = R(1e-6));

/// Integrates the functor @p integrand over the open interval (a, b) using the
/// double exponential rule.
/// @param integrand
///  The integrand that should be used for the integration.
/// @param a
///  The left endpoint of the integration interval.
/// @param b
///  The right endpoint of the integration interval.
/// @param maxIterations
///  The maximal number of iterations that the functions should use to achieve
///  the require precision. In each iteration the number of sampling points is
///  tripled.
/// @param epsilon
///  The upper bound on the integration error. If the function is not able to
///  achieved the required precision, an exception is thrown.
AURORA_API Real integrateDE(const std::function<Real (Real)>& integrand,
                            Real a, Real b, size_t maxIterations = 21,
                            Real epsilon = R(1e-6), Real hMax = R(4.3));

} // namespace Aurora

#endif
