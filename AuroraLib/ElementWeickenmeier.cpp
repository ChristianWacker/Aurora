//--- Aurora/AuroraLib/ElementWeickenmeier.cpp ---------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Element.hpp"
#include "AuroraLib/ElementHydrogen.hpp"
#include "AuroraLib/Exception.hpp"

namespace Aurora
{

/// Implements Element using the parameterization of Weickenmeier and Kohl
/// [Weickenmeier1991]
class AURORA_API ElementWeickenmeier : public Element
{
public:
  ElementWeickenmeier(size_t atomicNumber, int charge, Real v,
                      const std::array<Real, 6>& b);

  Real electronicFormFactor(Real kSqr, Real bFactor) const override;
  Real potential(Real rSqr, Real bFactor) const override;
  Real radiusSqr() const override;

private:
  std::array<Real, 6> mA;
  const std::array<Real, 6> mB;
};

ElementWeickenmeier::ElementWeickenmeier(size_t atomicNumber, int charge,
                                         Real v, const std::array<Real, 6>& b) :
  Element(atomicNumber, charge), mB(b)
{
  mA[0] = R(0.02395) / R(3.0) / (1 + v) * atomicNumber;
  mA[1] = mA[0];
  mA[2] = mA[0];
  mA[3] = v * mA[0];
  mA[4] = mA[3];
  mA[5] = mA[3];
}

Real ElementWeickenmeier::electronicFormFactor(Real kSqr, Real bFactor) const
{
  const Real sSqr = kSqr / (16 * Math::piSqr);

  Real result = 0;
  for (size_t i = 0; i < 6; ++i)
    result += mA[i] * (1 - std::exp(-mB[i] * sSqr));

  return result * std::exp(-bFactor * kSqr) / sSqr;
}

Real ElementWeickenmeier::potential(Real rSqr, Real bFactor) const
{
  const Real r = std::sqrt(rSqr);

  Real result = 0;
  for (size_t i = 0; i < 6; ++i)
  {
    const Real radicand = mB[i] + 16 * Math::piSqr * bFactor;
    result += mA[i] / r * (std::erf(R(0.5) * r / std::sqrt(bFactor)) -
                           std::erf(Math::twoPi * r / std::sqrt(radicand)));
  }

  return -8 * Math::piSqr * powerOf<2>(Phys::hBar) / Phys::mE * result;
}

Real ElementWeickenmeier::radiusSqr() const
{
  Real result = 0;
  for (size_t i = 0; i < 6; ++i)
    result += mA[i] * mB[i];

  return 3 * Phys::a0 * result / mAtomicNumber;
}

class AURORA_API FormFactorParametrizationWeickenmeier :
  public FormFactorParametrization
{
public:
  FormFactorParametrizationWeickenmeier()
  {
    // [Weickenmeier1991], table 1
    addElement(std::make_unique<ElementHydrogen>());
    addElement(ConstElementPtr(new ElementWeickenmeier( 2, 0, R(0.5), // HE
      {{R(2.542   ), R(8.743), R(1.269e1), R(4.371e-1), R(5.294   ), R(2.825e1)}})));
    addElement(ConstElementPtr(new ElementWeickenmeier( 6, 0, R(0.5), // C
      {{R(2.946e-1), R(3.934), R(2.498e1), R(2.528e1 ), R(2.547e1 ), R(4.670e1)}})));
    addElement(ConstElementPtr(new ElementWeickenmeier(14, 0, R(0.5), // SI
      {{R(1.737   ), R(3.043), R(3.057e1), R(5.070e-2), R(9.918e-1), R(8.618e1)}})));
    addElement(ConstElementPtr(new ElementWeickenmeier(92, 0, R(0.2), // U
      {{R(7.142e-2), R(1.149), R(9.212  ), R(9.592e-1), R(1.203   ), R(1.043e2)}})));
  }
};

class AURORA_API RegisterWeickenmeier
{
public:
  RegisterWeickenmeier()
  {
    FormFactorParametrization::
      registerClass<FormFactorParametrizationWeickenmeier>("Weickenmeier");
  }
} registerWeickenmeier;

} // namespace Aurora
