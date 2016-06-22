//--- Aurora/AuroraLib/Math.cpp ------------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MATH_HPP
#define AURORA_AURORA_LIB_MATH_HPP

#include "AuroraLib/Complex.hpp"
#include "AuroraLib/Exception.hpp"

#include <Eigen/Core>

namespace Aurora
{

/// mathematical constants in double precision
namespace MathDouble
{
  constexpr double pi       = 3.14159265358979323846;
  constexpr double twoPi    = 2.0 * pi;
  constexpr double fourPi   = 4.0 * pi;
  constexpr double halfPi   = 0.5 * pi;
  constexpr double thirdPi  = pi / 3.0;
  constexpr double piSqr    = pi * pi;
  constexpr double piCubed  = pi * pi * pi;
  constexpr double sqrt2    = 1.41421356237309504880;
  constexpr double degToRad = pi / 180.0;
}

/// mathematical constants
namespace Math
{
  constexpr Real pi       = static_cast<Real>(MathDouble::pi);
  constexpr Real twoPi    = static_cast<Real>(MathDouble::twoPi);
  constexpr Real fourPi   = static_cast<Real>(MathDouble::fourPi);
  constexpr Real halfPi   = static_cast<Real>(MathDouble::halfPi);
  constexpr Real thirdPi  = static_cast<Real>(MathDouble::thirdPi);
  constexpr Real piSqr    = static_cast<Real>(MathDouble::piSqr);
  constexpr Real piCubed  = static_cast<Real>(MathDouble::piCubed);
  constexpr Real sqrt2    = static_cast<Real>(MathDouble::sqrt2);
  constexpr Real degToRad = static_cast<Real>(MathDouble::degToRad);
}

typedef Eigen::Matrix<Real, 3, 3> Matrix3x3R;
typedef Eigen::Matrix<Real, 3, 1> Vector3R;

inline Vector3R min(const Vector3R& a, const Vector3R& b)
{
  return Vector3R(std::min(a.x(), b.x()),
                  std::min(a.y(), b.y()),
                  std::min(a.z(), b.z()));
}

inline Vector3R max(const Vector3R& a, const Vector3R& b)
{
  return Vector3R(std::max(a.x(), b.x()),
                  std::max(a.y(), b.y()),
                  std::max(a.z(), b.z()));
}

inline Vector3R operator/(const Vector3R& lhs, const Vector3R& rhs)
{
  return Vector3R(lhs.x() / rhs.x(),
                  lhs.y() / rhs.y(),
                  lhs.z() / rhs.z());
}

class AURORA_API EMath : public Exception
{
public:
  /// @copydoc EInvalidParameter::EInvalidParameter
  EMath(const std::string& message, const std::string& function,
        const std::string& filename, size_t line);
};

/// Linear interpolation (lerp) between @p a and @p b
inline Real lerp(Real a, Real b, Real frac)
{
  return a + (b - a) * frac;
}

template<class T>
constexpr inline bool isEven(T x)
{
  return (x % 2) == 0;
}

/// Returns true, if x is a power of two
template<class T>
constexpr inline bool isPowerOfTwo(T x)
{
  // This code checks first, if x is zero.
  // Then the code uses the fact that subtracting one from any power of two
  // changes the entire bit pattern starting with the leading one.
  return (x != 0) && ((x & (x - 1)) == 0);
}

namespace Internal
{
  /// general, odd
  template<class T, int64_t n, bool even>
  struct PowerOf
  {
    static_assert(n > 0, "PowerOf is only defined for positive exponents");
    constexpr static T eval(T x)
    {
      return x * PowerOf<T, n - 1, true>::eval(x);
    }
  };

  /// general, even
  template<class T, int64_t n>
  struct PowerOf<T, n, true>
  {
    static_assert(n > 0, "PowerOf is only defined for positive exponents");
    constexpr static T eval(T x)
    {
      return sqr(PowerOf<T, n / 2, isEven(n / 2)>::eval(x));
    }

  private:
    constexpr static T sqr(T x)
    {
      return x * x;
    }
  };

  template<class T>
  struct PowerOf<T, 1, false>
  {
    constexpr static T eval(T x)
    {
      return x;
    }
  };
}

/// Template to calculate the power x^n if n is an integer. For small n this
/// template is faster than a call to pow(x, n).
template<int64_t n>
constexpr inline float powerOf(float x)
{
  return Internal::PowerOf<float, n, isEven(n)>::eval(x);
}

/// Template to calculate the power x^n if n is an integer. For small n this
/// template is faster than a call to pow(x, n).
template<int64_t n>
constexpr inline double powerOf(double x)
{
  return Internal::PowerOf<double, n, isEven(n)>::eval(x);
}

AURORA_CONSTEXPR inline size_t factorial(size_t n)
{
  size_t result = 1;
  for (size_t i = 2; i <= n; ++i)
    result *= i;
  return result;
}

/// Returns the binomial coefficient
///  @f[
///   \begin{pmatrix}\alpha\\k\end{pmatrix} = \
///    \frac{\alpha(\alpha-1)(\alpha-2)\cdots(\alpha-k+1)}{k(k-1)(k-2)\cdots 1}
///  @f]
AURORA_CONSTEXPR inline Real binomial(Real alpha, int64_t k)
{
  if (k < 0)
    return 0;

  if (0 == k)
    return 1;

  Real result = 1.0;
  for (int64_t i = 1; i <= k; ++i)
  {
    result *= alpha / i;
    alpha -= 1;
  }

  return result;
}

/// Returns the modified Bessel function of the first kind of order zero. This
/// function is defined for all x.
AURORA_API Real besselI0(Real x);
/// Returns the modified Bessel function of the first kind of order one. This
/// function is defined for all x.
AURORA_API Real besselI1(Real x);
/// Returns the modified Bessel function of the second kind of order zero.
/// This function is only defined for positive x.
AURORA_API Real besselK0(Real x);
/// Returns the modified Bessel function of the second kind of order one.
/// This function is only defined for positive x.
AURORA_API Real besselK1(Real x);

} // namespace Aurora

#endif
