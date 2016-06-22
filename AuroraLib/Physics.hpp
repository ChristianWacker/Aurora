//--- Aurora/AuroraLib/Physics.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_PHYSICS_HPP
#define AURORA_AURORA_LIB_PHYSICS_HPP

#include "AuroraLib/Math.hpp"

#include <cmath>

namespace Aurora
{

/// physical constants in SI units. These must be double, as float can only
/// represent 10^-39..10^+39.
namespace SI
{
  /// Speed of light in vacuum in m / s (exact)
  constexpr double c = 299792458.0;
  /// elementary charge in C [CODATA2010]
  constexpr double e = 1.602176565e-19;
  /// reduced Planck constant in J s [CODATA2010]
  constexpr double hBar = 1.054571726e-34;
  /// electron mass in kg [CODATA2010]
  constexpr double mE = 9.10938291e-31;
  /// vacuum permeability in H / m (exact)
  constexpr double mu0 = 4e-7 * MathDouble::pi;

  /// vacuum permittivity in F / m (exact, derived)
  constexpr double eps0 = 1.0 / (mu0 * c * c);
  /// Bohr radius in m (derived)
  constexpr double a0 = MathDouble::fourPi * eps0 * hBar * hBar / mE / e / e;
}

/// physical constants: We are using a system of measurement with
///  [length] = Angstrom
///  [energy] = eV
///  [time]   = Angstrom / c
///  [charge] = e
/// Hence, we don't need c or e as these are one.
namespace Phys
{
  /// electron mass in eV
  constexpr Real mE   = static_cast<Real>(SI::mE / SI::e * SI::c * SI::c);
  /// reduced Planck constant in Angstrom
  constexpr Real hBar = static_cast<Real>(SI::hBar * 1e10 * SI::c / SI::e);
  /// vacuum permittivity in V * Angstrom
  constexpr Real eps0 = static_cast<Real>(SI::eps0 / SI::e * 1e-10);
  /// Bohr radius in Angstrom
  constexpr Real a0   = static_cast<Real>(SI::a0 * 1e10);
}

/// Returns the speed in fractions of the speed of light.
/// @details
///  The function uses the formula
///    @f[
///      v = c \sqrt{1-\gamma^-2}
///    @f]
inline Real lorentzFactorToSpeed(Real lorentzFactor)
{
  assert(lorentzFactor >= R(1.0));
  return std::sqrt(R(1.0) - R(1.0) / powerOf<2>(lorentzFactor));
}

} // namespace Aurora

#endif
