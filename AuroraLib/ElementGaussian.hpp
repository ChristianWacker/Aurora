//--- Aurora/AuroraLib/ElementGaussian.hpp -------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_ELEMENT_GAUSSIAN_HPP
#define AURORA_AURORA_LIB_ELEMENT_GAUSSIAN_HPP

#include "AuroraLib/Element.hpp"
#include "AuroraLib/Math.hpp"
#include "AuroraLib/Physics.hpp"

#include <cmath>

namespace Aurora
{

/// Implements Element using the parameterization of Doyle Turner [Doyle1968]
template<size_t n>
class ElementGaussian : public Element
{
  static_assert(n > 0, "n must be positive.");
public:
  ElementGaussian(size_t atomicNumber, int charge,
                  const std::array<Real, n>& a,
                  const std::array<Real, n>& b);

  Real phaseShift(Real rho, Real z, Real bFactor) const override;
  Real projectedPotential(Real rhoSqr, Real bFactor) const override;
  Real electronicFormFactor(Real kSqr, Real bFactor) const override;
  Real potential(Real rSqr, Real bFactor) const override;
  Real hybridPotential(Real kappaSqr, Real z, Real bFactor) const override;
  Real hybridPhaseShift(Real kappaSqr, Real z, Real bFactor) const override;
  Real radiusSqr() const override;

private:
  const std::array<Real, n> mA;
  const std::array<Real, n> mB;
};

template<size_t n>
ElementGaussian<n>::ElementGaussian(size_t atomicNumber, int charge,
                                    const std::array<Real, n>& a,
                                    const std::array<Real, n>& b) :
  Element(atomicNumber, charge), mA(a), mB(b)
{ }

template<size_t n>
Real ElementGaussian<n>::projectedPotential(Real rhoSqr, Real bFactor) const
{
  Real result = 0.0;

  for (size_t i = 0; i < n; ++i)
  {
    const Real b = mB[i] + 16 * Math::piSqr * bFactor;
    result += mA[i] / b * std::exp(-4 * Math::piSqr / b * rhoSqr);
  }

  return -8 * powerOf<2>(Math::pi * Phys::hBar) / Phys::mE * result;
}

template<size_t n>
Real ElementGaussian<n>::phaseShift(Real rho, Real z, Real bFactor) const
{
  using namespace Math;
  const Real rhoSqr = rho * rho;

  Real result = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const Real b = mB[i] + 16 * piSqr * bFactor;
    result += mA[i] / b * std::exp(-4 * piSqr * rhoSqr / b )
      * std::erf(twoPi * z / std::sqrt(b));
  }

  return 4 * piSqr * result;
}

template<size_t n>
Real ElementGaussian<n>::electronicFormFactor(Real kSqr, Real bFactor) const
{
  Real result = 0;

  for (size_t i = 0; i < n; ++i)
  {
    Real b = mB[i] / (16 * Math::piSqr) + bFactor;
    result += mA[i] * std::exp(-b * kSqr);
  }

  return result;
}

template<size_t n>
Real ElementGaussian<n>::potential(Real rSqr, Real bFactor) const
{
  using namespace Math;
  using namespace Phys;

  Real result = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const Real b = mB[i] + 16 * piSqr * bFactor;
    result += mA[i] * std::pow(b, R(-1.5)) * std::exp(-4 * piSqr / b * rSqr);
  }

  return -16 * result * std::pow(pi, R(2.5)) * powerOf<2>(hBar) / mE;
}

template<size_t n>
Real ElementGaussian<n>::hybridPotential(Real kappaSqr, Real z,
                                         Real bFactor) const
{
  using namespace Math;

  Real result = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const Real b = mB[i] + 16 * piSqr * bFactor;
    result += 4 * mA[i] * std::sqrt(powerOf<3>(pi) / b)
      * std::exp(-b * kappaSqr / (16 * piSqr) - 4 * piSqr * z * z / b);
  }

  return result;
}

template<size_t n>
Real ElementGaussian<n>::hybridPhaseShift(Real kappaSqr, Real z,
                                          Real bFactor) const
{
  using namespace Math;

  Real result = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const Real b = mB[i] + 16 * piSqr * bFactor;
    result += mA[i] * std::exp(-b * kappaSqr / (16 * piSqr))
      * std::erf(twoPi * z / std::sqrt(b));
  }

  return R(0.25) * result / pi;
}

template<size_t n>
Real ElementGaussian<n>::radiusSqr() const
{
  Real result = 0;

  for (size_t i = 0; i < n; ++i)
    result += mA[i];

  return 3 * Phys::a0 * result / mAtomicNumber;
}

} // namespace Aurora

#endif
