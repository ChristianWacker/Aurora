//--- Aurora/AuroraLib/Complex.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_COMPLEX_HPP
#define AURORA_AURORA_LIB_COMPLEX_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <cmath>
#include <complex>

namespace Aurora
{

typedef std::complex<Real> Complex;

inline Real absSqr(const Complex& z)
{
  return z.real() * z.real() + z.imag() * z.imag();
}

inline Complex fromArg(Real arg)
{
  return Complex(std::cos(arg), std::sin(arg));
}

/// Returns the real part of the complex number.
inline Real real(const Complex& z)
{
  return z.real();
}

/// Returns the imaginary part of the complex number.
inline Real imag(const Complex& z)
{
  return z.imag();
}

/// Returns the absolute value of the complex number. This is the distance to
/// the origin in the complex plane.
inline Real abs(const Complex& z)
{
  return std::hypot(z.real(), z.imag());
}

inline Real arg(const Complex& z)
{
  return std::atan2(z.imag(), z.real());
}

#if 0
/// Class to represent a complex number.
/// @brief Replacement for std::complex as this class has obviously performance
///  issues. Therefore, we provide direct element access and the default
///  constructor does not initialize the object.
class AURORA_API Complex
{
public:
  /// Default constructor does not initialize the members.
  Complex() { }

  /// Initializes the complex number with a real number. Hence the imaginary
  /// part ist set to zero.
  constexpr explicit Complex(Real r) :
    r(r), i(Real())
  { }

  /// Initializes the complex number with its real and imaginary part.
  constexpr Complex(Real r, Real i) :
    r(r), i(i)
  { }

  explicit Complex(std::complex<Real> c) :
    r(c.real()), i(c.imag())
  { }

  /// Assigns a new real number. The imaginary part is set to zero.
  Complex& operator=(Real rhs)
  {
    r = rhs;
    i = Real();
    return *this;
  }

  /// Adds to this complex number another complex number and returns a
  /// reference to the result.
  Complex& operator+=(const Complex& rhs)
  {
    r += rhs.r;
    i += rhs.i;
    return *this;
  }

  /// Subtracts form this complex number another complex number and returns a
  /// reference to the result.
  Complex& operator-=(const Complex& rhs)
  {
    r -= rhs.r;
    i -= rhs.i;
    return *this;
  }

  /// Multiplies this complex number with another complex number and returns a
  /// reference to the result.
  Complex& operator*=(const Complex& rhs)
  {
    const Real savedR = r;
    r = r * rhs.r - i * rhs.i;
    i = savedR * rhs.i + i * rhs.r;
    return *this;
  }

  /// Multiplies this complex number with a real number and returns a
  /// reference to the result.
  Complex& operator*=(Real rhs)
  {
    r *= rhs;
    i *= rhs;
    return *this;
  }

  /// Divides this complex number by another complex number and returns a
  /// reference to the result.
  Complex& operator/=(const Complex& rhs)
  {
    const Real denominator = rhs.r * rhs.r + rhs.i * rhs.i;
    const Real savedR = r;
    r = (r * rhs.r + i * rhs.i) / denominator;
    i = (i * rhs.r - savedR * rhs.i) / denominator;
    return *this;
  }

  /// Divides this complex number by a real number and returns a reference to
  /// the result.
  Complex& operator/=(Real rhs)
  {
    r /= rhs;
    i /= rhs;
    return *this;
  }

  /// @name Unary operators
  /// @{
  /// Returns this complex number.
  Complex operator+() const
  {
    return *this;
  }

  /// Returns this complex number multiplied by -1.
  Complex operator-() const
  {
    return Complex(-r, -i);
  }
  /// @}

  /// @name Arithmetic operators
  /// @{
  /// Returns the sum of a real and a complex number.
  friend Complex operator+(Real lhs, const Complex& rhs)
  {
    return Complex(lhs + rhs.r, rhs.i);
  }

  /// Returns the sum of a complex and a real number.
  Complex operator+(Real rhs) const
  {
    return Complex(r + rhs, i);
  }

  /// Returns the difference of real number and a complex number.
  friend Complex operator-(Real lhs, const Complex& rhs)
  {
    return Complex(lhs - rhs.r, -rhs.i);
  }

  /// Returns the difference of a complex number and a real number.
  Complex operator-(Real rhs) const
  {
    return Complex(r - rhs, i);
  }

  /// Returns the product of a real number and a complex number.
  friend Complex operator*(Real lhs, const Complex& rhs)
  {
    return Complex(lhs * rhs.r, lhs * rhs.i);
  }

  /// Returns the product of a complex number and a real number.
  Complex operator*(Real rhs) const
  {
    return Complex(r * rhs, i * rhs);
  }

  /// Returns the quotient of a real number and a complex number.
  friend Complex operator/(Real lhs, const Complex& rhs)
  {
    const Real denominator = rhs.r * rhs.r + rhs.i * rhs.i;
    return Complex((lhs * rhs.r) / denominator,
                   (-lhs * rhs.i) / denominator);
  }

  /// Returns the quotient of a complex number and a real number.
  Complex operator/(Real rhs) const
  {
    return Complex(r / rhs, i / rhs);
  }
  /// @}

  /// Creates a new unit complex number with the given argument. This call is
  /// equivalent to @ref fromPolar(1.0, arg).
  static Complex fromArg(Real arg)
  {
    return Complex(std::cos(arg), std::sin(arg));
  }

  /// Returns a new complex form its representation in polar coordinates.
  static Complex fromPolar(Real abs, Real arg)
  {
    return Complex(abs * std::cos(arg), abs * std::sin(arg));
  }

  /// Conversion
  explicit operator std::complex<Real>() const
  {
    return std::complex<Real>(r, i);
  }

  Real r, i;
};

/// Returns the sum of two complex numbers.
inline Complex operator+(const Complex& lhs, const Complex& rhs)
{
  return Complex(lhs.r + rhs.r, lhs.i + rhs.i);
}

/// Returns the difference of two complex numbers.
inline Complex operator-(const Complex& lhs, const Complex& rhs)
{
  return Complex(lhs.r - rhs.r, lhs.i - rhs.i);
}

/// Returns the product of two complex numbers.
inline Complex operator*(const Complex& lhs, const Complex& rhs)
{
  return Complex(lhs.r * rhs.r - lhs.i * rhs.i, lhs.r * rhs.i + lhs.i * rhs.r);
}

/// Returns the quotient of two complex numbers.
inline Complex operator/(const Complex& lhs, const Complex& rhs)
{
  const Real denominator = rhs.r * rhs.r + rhs.i * rhs.i;
  return Complex((lhs.r * rhs.r + lhs.i * rhs.i) / denominator,
                 (lhs.i * rhs.r - lhs.r * rhs.i) / denominator);
}

/// Checks two complex number for equality.
inline bool operator==(const Complex& lhs, const Complex& rhs)
{
  return (lhs.r == rhs.r) && (lhs.i == rhs.i);
}

/// Checks two complex number for inequality.
inline bool operator!=(const Complex& lhs, const Complex& rhs)
{
  return (lhs.r != rhs.r) || (lhs.i != rhs.i);
}

/// Returns the complex conjugated number.
inline Complex conj(const Complex& z)
{
  return Complex(z.r, -z.i);
}

/// Returns the real part of the complex number.
inline Real real(const Complex& z)
{
  return z.r;
}

/// Returns the imaginary part of the complex number.
inline Real imag(const Complex& z)
{
  return z.i;
}

/// Returns the absolute value of the complex number. This is the distance to
/// to the origin in the complex plane.
inline Real abs(const Complex& z)
{
  return std::hypot(z.r, z.i);
}

/// Returns the squared absolute value of the complex number.
inline Real absSqr(const Complex& z)
{
  return z.r * z.r + z.i * z.i;
}

/// Returns the argument of the complex number. This is the angle in complex
/// plane with positive x axis. The return value is in the interval [-pi, pi].
inline Real arg(const Complex& z)
{
  return std::atan2(z.i, z.r);
}

/// Streams the complex number into an std::ostream.
inline std::ostream& operator<<(std::ostream& stream, const Complex& c)
{
  stream << '(' << c.r << ", " << c.i << ')';
  return stream;
}

/// Suffix for imaginary numbers.
constexpr Complex operator"" _i(long double imag)
{
  return Complex(0.0, static_cast<Real>(imag));
}
#endif

} // namespace Aurora

#endif
