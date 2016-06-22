//--- Aurora/AuroraLib/ElementHydrogen.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_ELEMENT_HYDROGEN_HPP
#define AURORA_AURORA_LIB_ELEMENT_HYDROGEN_HPP

#include "AuroraLib/Element.hpp"
#include "AuroraLib/Math.hpp"
#include "AuroraLib/Physics.hpp"

namespace Aurora
{

/// The exact form factor of the hydrogen atom
class AURORA_API ElementHydrogen : public Element
{
public:
  ElementHydrogen();

  Real projectedPotential(Real rhoSqr, Real bFactor) const override;
  Real electronicFormFactor(Real kSqr, Real bFactor) const override;
  Real potential(Real rSqr, Real bFactor) const override;
  Real radiusSqr() const override;
};

/// The exact form factor of the hydrogen ion
class AURORA_API ElementHydrogenIon : public Element
{
public:
  ElementHydrogenIon();

  Real phaseShift(Real rho, Real z, Real bFactor) const override;
  Real electronicFormFactor(Real kSqr, Real bFactor) const override;
  Real potential(Real rSqr, Real bFactor) const override;
};

inline ElementHydrogen::ElementHydrogen() :
  Element(1, 0)
{}

inline Real ElementHydrogen::projectedPotential(Real rhoSqr, Real) const
{
  const Real rho = std::sqrt(rhoSqr);
  const Real x = 2 * rho / Phys::a0;
  Real result = besselK0(x) + besselK1(x) / x;
  return -powerOf<2>(Phys::hBar) / (Phys::mE * Phys::a0) * result;
}

inline Real ElementHydrogen::electronicFormFactor(Real kSqr, Real bFactor) const
{
  constexpr Real a0Sqr = powerOf<2>(Phys::a0);
  return (R(2.0) - R(32.0) / powerOf<2>(R(4.0) + a0Sqr * kSqr))
         * std::exp(-bFactor * kSqr) / (Phys::a0 * kSqr);
}

inline Real ElementHydrogen::potential(Real rSqr, Real bFactor) const
{
  (void)bFactor;
  assert(0 == bFactor && "Non-zero bFactor not implemented.");

  using namespace Phys;

  const Real r = std::sqrt(rSqr);
  return -hBar * hBar * (a0 + r) * std::exp(-2 * r / a0) / (mE * r * powerOf<2>(a0));
}

inline Real ElementHydrogen::radiusSqr() const
{
  return 3 * powerOf<2>(Phys::a0) / mAtomicNumber;
}

inline ElementHydrogenIon::ElementHydrogenIon() :
  Element(1, 1)
{}

inline Real ElementHydrogenIon::electronicFormFactor(Real, Real) const
{
  return std::numeric_limits<Real>::quiet_NaN();
}

inline Real ElementHydrogenIon::phaseShift(Real rho, Real z, Real) const
{
  return std::asinh(z / rho) / Phys::a0;
}

inline Real ElementHydrogenIon::potential(Real rSqr, Real bFactor) const
{
  (void)bFactor;
  assert(0 == bFactor && "Non-zero bFactor not implemented.");

  const Real r = std::sqrt(rSqr);
  return -R(1.0) / (Math::fourPi * Phys::eps0 * r);
}

} // namespace Aurora

#endif
