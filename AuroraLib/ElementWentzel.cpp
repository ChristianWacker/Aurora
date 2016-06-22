//--- Aurora/AuroraLib/ElementWentzel.cpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Element.hpp"
#include "AuroraLib/Exception.hpp"
#include "AuroraLib/Physics.hpp"

namespace Aurora
{

/// Implements Elements using the Wentzel(Yukawa) potential
class AURORA_API ElementWentzel : public Element
{
public:
  ElementWentzel(size_t atomicNumber);

  Real electronicFormFactor(Real kSqr, Real bFactor) const override;
  Real potential(Real rSqr, Real bFactor) const override;
  Real xRayFormFactor(Real kSqr) const override;
  Real radiusSqr() const override;

private:
  // the screening radius
  const Real mR0;
};

class AURORA_API FormFactorParametrizationWentzel :
  public FormFactorParametrization
{
public:
  FormFactorParametrizationWentzel()
  {
    for (size_t i = 1; i <= Element::mNumElements; ++i)
      addElement(std::make_unique<ElementWentzel>(i));
  }
};

ElementWentzel::ElementWentzel(size_t atomicNumber) :
  Element(atomicNumber, 0),
  mR0(Phys::a0 * std::pow(static_cast<Real>(atomicNumber), -R(1.0) / R(3.0)))
{}

Real ElementWentzel::electronicFormFactor(Real kSqr, Real bFactor) const
{
  return 2 * atomicNumber() / (Phys::a0 * (kSqr + 1 / (mR0 * mR0)))
    * std::exp(-bFactor * kSqr);
}

Real ElementWentzel::potential(Real rSqr, Real bFactor) const
{
  using namespace Phys;

  const Real r = std::sqrt(rSqr);
  const Real sqrtB = std::sqrt(bFactor);

  if (0 == bFactor)
    return -hBar * hBar * atomicNumber() * std::exp(-r / mR0) / (a0 * mE * r);

  Real result =
      std::exp(-r / mR0) * std::erfc(sqrtB / mR0 - R(0.5) * r / sqrtB)
    - std::exp( r / mR0) * std::erfc(sqrtB / mR0 + R(0.5) * r / sqrtB);

  result *= -hBar * hBar * atomicNumber() * std::exp(bFactor / powerOf<2>(mR0));
  result /= 2 * mE * a0 * r;

  return result;
}

Real ElementWentzel::xRayFormFactor(Real kSqr) const
{
  return mAtomicNumber / (1 + kSqr * mR0 * mR0);
}

Real ElementWentzel::radiusSqr() const
{
  return 6 * mR0 * mR0;
}

class AURORA_API RegisterWentzel
{
public:
  RegisterWentzel()
  {
    FormFactorParametrization::registerClass<FormFactorParametrizationWentzel>
      ("Wentzel");
  }
} registerWentzel;

} // namespace Aurora
