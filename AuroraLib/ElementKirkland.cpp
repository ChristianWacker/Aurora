//--- Aurora/AuroraLib/ElementKirkland.cpp -------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Element.hpp"
#include "AuroraLib/Exception.hpp"
#include "AuroraLib/Math.hpp"
#include "AuroraLib/Physics.hpp"

namespace Aurora
{

/// Implements Element using the parameterization of Kirkland [Kirkland2010]
class AURORA_API ElementKirkland : public Element
{
public:
  // number of Lorentz functions
  static const size_t numLorentzian = 3;
  // number of Gauss functions
  static const size_t numGaussian = 3;

  ElementKirkland(size_t atomicNumber, int charge,
                  const std::array<Real, numLorentzian>& a,
                  const std::array<Real, numLorentzian>& b,
                  const std::array<Real, numGaussian>& c,
                  const std::array<Real, numGaussian>& d);

  Real projectedPotential(Real rhoSqr, Real bFactor) const override;
  Real electronicFormFactor(Real kSqr, Real bFactor) const override;
  Real potential(Real rSqr, Real bFactor) const override;
  Real hybridPotential(Real kappaSqr, Real z, Real bFactor) const override;
  Real hybridPhaseShift(Real kappaSqr, Real z, Real bFactor) const override;
  Real radiusSqr() const override;

private:
  const std::array<Real, numLorentzian> mA;
  const std::array<Real, numLorentzian> mB;
  const std::array<Real, numGaussian> mC;
  const std::array<Real, numGaussian> mD;
};

class AURORA_API FormFactorParametrizationKirkland :
  public FormFactorParametrization
{
public:
  FormFactorParametrizationKirkland();
};

ElementKirkland::ElementKirkland(size_t atomicNumber, int charge,
                                 const std::array<Real, numLorentzian>& a,
                                 const std::array<Real, numLorentzian>& b,
                                 const std::array<Real, numGaussian>& c,
                                 const std::array<Real, numGaussian>& d) :
  Element(atomicNumber, charge), mA(a), mB(b), mC(c), mD(d)
{}

Real ElementKirkland::projectedPotential(Real rhoSqr, Real bFactor) const
{
  (void)bFactor;
  assert(0 == bFactor && "b-Factor is ignored.");
  // todo: bFactor
  const Real rho = std::sqrt(rhoSqr);

  Real result = 0;
  for (size_t i = 0; i < numLorentzian; ++i)
    result += 2 * mA[i] * besselK0(Math::twoPi * rho * std::sqrt(mB[i]));

  for (size_t i = 0; i < numGaussian; ++i)
    result += mC[i] / mD[i] * std::exp(-Math::piSqr * rhoSqr / mD[i]);

  return -2 * Math::piSqr * powerOf<2>(Phys::hBar) / Phys::mE * result;
}

Real ElementKirkland::electronicFormFactor(Real kSqr, Real bFactor) const
{
  const Real x = kSqr / (4 * Math::piSqr);

  Real result = 0;
  for (size_t i = 0; i < numLorentzian; ++i)
    result += mA[i] / (mB[i] + x);

  for (size_t i = 0; i < numGaussian; ++i)
    result += mC[i] * std::exp(-mD[i] * x);

  return result * std::exp(-bFactor * kSqr);
}

Real ElementKirkland::potential(Real rSqr, Real bFactor) const
{
  using namespace Math;
  using namespace Phys;

  const Real r = std::sqrt(rSqr);
  Real result = 0;

  if (0 == bFactor)
  {
    // We can use simpler expressions
    for (size_t i = 0; i < numLorentzian; ++i)
      result += piSqr * mA[i] / r * std::exp(-twoPi * std::sqrt(mB[i]) * r);

    for (size_t i = 0; i < numGaussian; ++i)
    {
      result += std::pow(pi, R(2.5)) * mC[i] * std::pow(mD[i], R(-1.5)) *
        std::exp(-piSqr * rSqr / mD[i]);
    }

    return -2 * powerOf<2>(hBar) / mE * result;
  }

  const Real sqrtBFactor2 = 2 * std::sqrt(bFactor);
  for (size_t i = 0; i < numLorentzian; ++i)
  {
    result -= mA[i] * powerOf<2>(hBar) * piSqr / (mE * r)
     * std::exp(4 * piSqr * mB[i] * bFactor)
     * (std::exp(-twoPi * std::sqrt(mB[i]) * r) *
        std::erfc(twoPi * std::sqrt(mB[i] * bFactor) - r / sqrtBFactor2)
        - std::exp(twoPi * std::sqrt(mB[i]) * r) *
        std::erfc(twoPi * std::sqrt(mB[i] * bFactor) + r / sqrtBFactor2));
  }

  for (size_t i = 0; i < numGaussian; ++i)
  {
    const Real x = mD[i] + 4 * piSqr * bFactor;
    result -= 2 * powerOf<2>(hBar) / mE * std::pow(pi, R(2.5))
     * mC[i] * std::pow(x, -R(1.5)) * std::exp(-piSqr * rSqr / x);
  }

  return result;
}

Real ElementKirkland::hybridPotential(Real kappaSqr, Real z, Real bFactor) const
{
  using namespace Math;
  using namespace Phys;

  Real result = 0;

  if (0 == bFactor)
  {
    for (size_t i = 0; i < numLorentzian; ++i)
    {
      Real a = std::sqrt(kappaSqr + 4* piSqr * mB[i]);
      result += 4 * powerOf<3>(pi) * mA[i] / a * std::exp(-a * abs(z));
    }
  }
  else
  {
    const Real sqrtBFactor = std::sqrt(bFactor);

    for (size_t i = 0; i < numLorentzian; ++i)
    {
      const Real aSqr = kappaSqr + 4 * piSqr * mB[i];
      const Real a = std::sqrt(aSqr);

      const Real g1 = std::erfc(a * sqrtBFactor + R(0.5) * z / sqrtBFactor);
      const Real g2 = std::erfc(a * sqrtBFactor - R(0.5) * z / sqrtBFactor);

      Real gsum = 0;
      if (0 != g1)
        gsum += std::exp( a * z) * g1;

      if (0 != g2)
        gsum += std::exp(-a * z) * g2;

      result += 2 * piCubed * mA[i] * std::exp((aSqr - kappaSqr) * bFactor)
        / a * gsum;
    }
  }

  for (size_t i = 0; i < numGaussian; ++i)
  {
    const Real d = mD[i] + 4 * piSqr * bFactor;
    result += 2 * mC[i] * std::sqrt(piCubed / d)
     * std::exp(-d * kappaSqr / (4 * piSqr) - piSqr * z * z / d);
  }

  return result;
}

Real ElementKirkland::hybridPhaseShift(Real kappaSqr, Real z,
                                       Real bFactor) const
{
  using namespace Math;

  Real result = 0;
  if (0 == bFactor)
  {
    for (size_t i = 0; i < numLorentzian; ++i)
    {
      const Real aSqr = kappaSqr + 4 * piSqr * mB[i];
      const Real a = std::sqrt(aSqr);
      result += pi * mA[i] / aSqr * std::copysign(1 - std::exp(-a * abs(z)), z);
    }
  }
  else
  {
    const Real sqrtBFactor = std::sqrt(bFactor);

    for (size_t i = 0; i < numLorentzian; ++i)
    {
      const Real aSqr = kappaSqr + 4 * piSqr * mB[i];
      const Real a = std::sqrt(aSqr);

      const Real g1 = std::erfc(a * sqrtBFactor + R(0.5) * z / sqrtBFactor);
      const Real g2 = std::erfc(a * sqrtBFactor - R(0.5) * z / sqrtBFactor);

      // To avoid NANs by multiplying 0 times inf, we simply ignore the
      // problematic term, if it has been already evaluated to be zero.
      Real gsum = 0;
      if (0 != g1)
        gsum += std::exp( a * z) * g1;

      if (0 != g2)
        gsum -= std::exp(-a * z) * g2;

      result += pi * mA[i] * std::exp(-bFactor * kappaSqr) / aSqr
        * (R(1.0) - std::erfc(R(0.5) * z / sqrtBFactor)
           + R(0.5) * std::exp(aSqr * bFactor) * gsum);
    }
  }

  for (size_t i = 0; i < numGaussian; ++i)
  {
    Real d = mD[i] + 4 * piSqr * bFactor;
    result += R(0.25) * mC[i] / pi * std::exp(-R(0.25) * d * kappaSqr / piSqr)
     * std::erf(pi * z / std::sqrt(d));
  }

  return result;
}

Real ElementKirkland::radiusSqr() const
{
  Real result = 0;

  for (size_t i = 0; i < numLorentzian; ++i)
    result += mA[i] / mB[i];

  for (size_t i = 0; i < numGaussian; ++i)
    result += mC[i];

  return result;
}

FormFactorParametrizationKirkland::FormFactorParametrizationKirkland()
{
  // [Kirkland2010]), R(p. 253 ff
  addElement(ConstElementPtr(new ElementKirkland(1, 0, // H
    {{R(4.20298324e-3), R(6.27762505e-2), R(3.00907347e-2)}},
    {{R(2.25350888e-1), R(2.25366950e-1), R(2.25331756e-1)}},
    {{R(6.77756695e-2), R(3.56609237e-3), R(2.76135815e-2)}},
    {{R(4.38854001e-0), R(4.03884823e-1), R(1.44490166e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(2, 0, // HE
    {{R(1.87543704e-5), R(4.10595800e-4), R(1.96300059e-1)}},
    {{R(2.12427997e-1), R(3.32212279e-1), R(5.17325152e-1)}},
    {{R(8.36015738e-3), R(2.95102022e-2), R(4.65928982e-7)}},
    {{R(3.66668239e-1), R(1.37171827e-0), R(3.75768025e-4)}})));
  addElement(ConstElementPtr(new ElementKirkland(3, 0, // LI
    {{R(7.45843816e-2), R(7.15382250e-2), R(1.45315229e-1)}},
    {{R(8.81151424e-1), R(4.59142904e-2), R(8.81301714e-1)}},
    {{R(1.12125769e-0), R(2.51736525e-3), R(3.58434971e-1)}},
    {{R(1.88483665e-1), R(1.59189995e-1), R(6.12371000e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(4, 0, // BE
    {{R(6.11642897e-2), R(1.25755034e-1), R(2.00831548e-1)}},
    {{R(9.90182132e-2), R(9.90272412e-2), R(1.87392509e-0)}},
    {{R(7.87242876e-1), R(1.58847850e-3), R(2.73962031e-1)}},
    {{R(9.32794929e-0), R(8.91900236e-2), R(3.20687658e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(5, 0, // B
    {{R(1.25716066e-1), R(1.73314452e-1), R(1.84774811e-1)}},
    {{R(1.48258830e-1), R(1.48257216e-1), R(3.34227311e-0)}},
    {{R(1.95250221e-1), R(5.29642075e-1), R(1.08230500e-3)}},
    {{R(1.97339463e-0), R(5.70035553e-0), R(5.64857237e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(6, 0, // C
    {{R(2.12080767e-1), R(1.99811865e-1), R(1.68254385e-1)}},
    {{R(2.08605417e-1), R(2.08610186e-1), R(5.57870773e-0)}},
    {{R(1.42048360e-1), R(3.63830672e-1), R(8.35012044e-4)}},
    {{R(1.33311887e-0), R(3.80800263e-0), R(4.03982620e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(7, 0, // N
    {{R(5.33015554e-1), R(5.29008883e-2), R(9.24159648e-2)}},
    {{R(2.90952515e-1), R(1.03547896e-1), R(1.03540028e-1)}},
    {{R(2.61799101e-1), R(8.80262108e-4), R(1.10166555e-1)}},
    {{R(2.76252723e-0), R(3.47681236e-2), R(9.93421736e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(8, 0, // O
    {{R(3.39969204e-1), R(3.07570172e-1), R(1.30369072e-1)}},
    {{R(3.81570280e-1), R(3.81571436e-1), R(1.91919745e-1)}},
    {{R(8.83326058e-2), R(1.96586700e-1), R(9.96220028e-4)}},
    {{R(7.60635525e-1), R(2.07401094e-0), R(3.03266869e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(9, 0, // F
    {{R(2.30560593e-1), R(5.26889648e-1), R(1.24346755e-1)}},
    {{R(4.80754213e-1), R(4.80763895e-1), R(3.95306720e-1)}},
    {{R(1.24616894e-3), R(7.20452555e-2), R(1.53075777e-1)}},
    {{R(2.62181803e-2), R(5.92495593e-1), R(1.59127671e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(10, 0, // NE
    {{R(4.08371771e-1), R(4.54418858e-1), R(1.44564923e-1)}},
    {{R(5.88228627e-1), R(5.88288655e-1), R(1.21246013e-2)}},
    {{R(5.91531395e-2), R(1.24003718e-1), R(1.64986037e-3)}},
    {{R(4.63963540e-1), R(1.23413025e-0), R(2.05869217e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(11, 0, // NA
    {{R(1.36471662e-1), R(7.70677865e-1), R(1.56862014e-1)}},
    {{R(4.99965301e-2), R(8.81899664e-1), R(1.61768579e-1)}},
    {{R(9.96821513e-1), R(3.80304670e-2), R(1.27685089e-1)}},
    {{R(2.00132610e-1), R(2.60516254e-1), R(6.99559329e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(12, 0, // MG
    {{R(3.04384121e-1), R(7.56270563e-1), R(1.01164809e-1)}},
    {{R(8.42014377e-2), R(1.64065598e-0), R(2.97142975e-1)}},
    {{R(3.45203403e-2), R(9.71751327e-1), R(1.20593012e-1)}},
    {{R(2.16596094e-1), R(1.21236852e-1), R(5.60865838e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(13, 0, // AL
    {{R(7.77419424e-1), R(5.78312036e-2), R(4.26386499e-1)}},
    {{R(2.71058227e-0), R(7.17532098e-1), R(9.13331555e-2)}},
    {{R(1.13407220e-1), R(7.90114035e-1), R(3.23293496e-2)}},
    {{R(4.48867451e-1), R(8.66366718e-0), R(1.78503463e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(14, 0, // SI
    {{R(1.06543892e-0), R(1.20143691e-1), R(1.80915263e-1)}},
    {{R(1.04118455e-0), R(6.87113368e-1), R(8.87533926e-2)}},
    {{R(1.12065620e-0), R(3.05452816e-2), R(1.59963502e-0)}},
    {{R(3.70062619e-0), R(2.14097897e-1), R(9.99096638e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(15, 0, // P
    {{R(1.05284447e-0), R(2.99440284e-1), R(1.17460748e-1)}},
    {{R(1.31962590e-0), R(1.28460520e-1), R(1.02190163e-2)}},
    {{R(9.60643452e-1), R(2.63555748e-2), R(1.38059330e-0)}},
    {{R(2.87477555e-0), R(1.82076844e-1), R(7.49165526e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(16, 0, // S
    {{R(1.01646916e-0), R(4.41766748e-1), R(1.21503863e-1)}},
    {{R(1.69181965e-0), R(1.74180288e-1), R(1.67011091e-2)}},
    {{R(8.27966670e-1), R(2.33022533e-2), R(1.18302846e-0)}},
    {{R(2.30342810e-0), R(1.56954150e-1), R(5.85782891e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(17, 0, // CL
    {{R(9.44221116e-1), R(4.37322049e-1), R(2.54547926e-1)}},
    {{R(2.40052374e-1), R(9.30510439e-0), R(9.30486346e-0)}},
    {{R(5.47763323e-2), R(8.00087488e-1), R(1.07488641e-2)}},
    {{R(1.68655688e-1), R(2.97849774e-0), R(6.84240646e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(18, 0, // AR
    {{R(1.06983288e-0), R(4.24631786e-1), R(2.43897949e-1)}},
    {{R(2.87791022e-1), R(1.24156957e-1), R(1.24158868e-1)}},
    {{R(4.79446296e-2), R(7.64958952e-1), R(8.23128431e-3)}},
    {{R(1.36979796e-1), R(2.43940729e-0), R(5.27258749e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(19, 0, // K
    {{R(6.92717865e-1), R(9.65161085e-1), R(1.48466588e-1)}},
    {{R(7.10849990e-0), R(3.57532901e-1), R(3.93763275e-2)}},
    {{R(2.64645027e-2), R(1.80883768e-0), R(5.43900018e-1)}},
    {{R(1.03591321e-1), R(3.22845199e-1), R(1.67791374e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(20, 0, // CA
    {{R(3.66902871e-1), R(8.66378999e-1), R(6.67203300e-1)}},
    {{R(6.14274129e-2), R(5.70881727e-1), R(7.82965639e-0)}},
    {{R(4.87743636e-1), R(1.82406314e-0), R(2.20248453e-2)}},
    {{R(1.32531318e-0), R(2.10056032e-1), R(9.11853450e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(21, 0, // SC
    {{R(3.78871777e-1), R(9.00022505e-1), R(7.15288914e-1)}},
    {{R(6.98910162e-2), R(5.21061541e-1), R(7.87707920e-0)}},
    {{R(1.88640973e-2), R(4.07945949e-1), R(1.61786540e-0)}},
    {{R(8.17512708e-2), R(1.11141388e-0), R(1.80840759e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(22, 0, // TI
    {{R(3.62383267e-1), R(9.84232966e-1), R(7.41715642e-1)}},
    {{R(7.54707114e-2), R(4.97757309e-1), R(8.17659391e-0)}},
    {{R(3.62555269e-1), R(1.49159390e-0), R(1.61659509e-2)}},
    {{R(9.55524906e-1), R(1.62221677e-1), R(7.33140839e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(23, 0, // V
    {{R(3.52961378e-1), R(7.46791014e-1), R(1.08364068e-0)}},
    {{R(8.19204103e-2), R(8.81189511e-0), R(5.10646075e-1)}},
    {{R(1.39013610e-0), R(3.31273356e-1), R(1.40422612e-2)}},
    {{R(1.48901841e-1), R(8.38543079e-1), R(6.57432678e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(24, 0, // CR
    {{R(1.34348379e-0), R(5.07040328e-1), R(4.26358955e-1)}},
    {{R(1.25814353e-0), R(1.15042811e-1), R(8.53660389e-2)}},
    {{R(1.17241826e-2), R(5.11966516e-1), R(3.38285828e-1)}},
    {{R(6.00177061e-2), R(1.53772451e-0), R(6.62418319e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(25, 0, // MN
    {{R(3.26697613e-1), R(7.17297000e-1), R(1.33212464e-0)}},
    {{R(8.88813083e-2), R(1.11300198e-1), R(5.82141104e-1)}},
    {{R(2.80801702e-1), R(1.15499241e-0), R(1.11984488e-2)}},
    {{R(6.71583145e-1), R(1.26825395e-1), R(5.32334467e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(26, 0, // FE
    {{R(3.13454847e-1), R(6.89290016e-1), R(1.47141531e-0)}},
    {{R(8.99325756e-2), R(1.30366038e-1), R(6.33345291e-1)}},
    {{R(1.03298688e-0), R(2.58280285e-1), R(1.03460690e-2)}},
    {{R(1.16783425e-1), R(6.09116446e-1), R(4.81610627e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(27, 0, // CO
    {{R(3.15878278e-1), R(1.60139005e-0), R(6.56394338e-1)}},
    {{R(9.46683246e-2), R(6.99436449e-1), R(1.56954403e-1)}},
    {{R(9.36746624e-1), R(9.77562646e-3), R(2.38378578e-1)}},
    {{R(1.09392410e-1), R(4.37446816e-2), R(5.56286483e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(28, 0, // NI
    {{R(1.72254630e-0), R(3.29543044e-1), R(6.23007200e-1)}},
    {{R(7.76606908e-1), R(1.02262360e-1), R(1.94156207e-1)}},
    {{R(9.43496513e-3), R(8.54063515e-1), R(2.21073515e-1)}},
    {{R(3.98684596e-2), R(1.04078166e-1), R(5.10869330e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(29, 0, // CU
    {{R(3.58774531e-1), R(1.76181348e-0), R(6.36905053e-1)}},
    {{R(1.06153463e-1), R(1.01640995e-0), R(1.53659093e-1)}},
    {{R(7.44930667e-3), R(1.89002347e-1), R(2.29619589e-1)}},
    {{R(3.85345989e-2), R(3.98427790e-1), R(9.01419843e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(30, 0, // ZN
    {{R(5.70893973e-1), R(1.98908856e-0), R(3.06060585e-1)}},
    {{R(1.26534614e-1), R(2.17781965e-0), R(3.78619003e-1)}},
    {{R(2.35600223e-1), R(3.97061102e-1), R(6.85657228e-3)}},
    {{R(3.67019041e-1), R(8.66419596e-1), R(3.35778823e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(31, 0, // GA
    {{R(6.25528464e-1), R(2.05302901e-0), R(2.89608120e-1)}},
    {{R(1.10005650e-1), R(2.41095786e-0), R(4.78685736e-1)}},
    {{R(2.07910594e-1), R(3.45079617e-1), R(6.55634298e-3)}},
    {{R(3.27807224e-1), R(7.43139061e-1), R(3.09411369e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(32, 0, // GE
    {{R(5.90952690e-1), R(5.39980660e-1), R(2.00626188e-0)}},
    {{R(1.18375976e-1), R(7.18937433e-1), R(1.39304889e-0)}},
    {{R(7.49705041e-1), R(1.83581347e-1), R(9.52190743e-3)}},
    {{R(6.89943350e-0), R(3.64667232e-1), R(2.69888650e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(33, 0, // AS
    {{R(7.77875218e-1), R(5.93848150e-1), R(1.95918751e-0)}},
    {{R(1.50733157e-1), R(1.42882209e-2), R(1.74750339e-0)}},
    {{R(1.79880226e-1), R(8.63267222e-1), R(9.59053427e-3)}},
    {{R(3.31800852e-1), R(5.85490274e-0), R(2.33777569e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(34, 0, // SE
    {{R(9.58390681e-1), R(6.03851342e-1), R(1.90828931e-0)}},
    {{R(1.83775557e-1), R(1.96819224e-2), R(2.15082053e-0)}},
    {{R(1.73885956e-1), R(9.35265145e-1), R(8.62254658e-3)}},
    {{R(3.00006024e-1), R(4.92471215e-0), R(2.12308108e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(35, 0, // BR
    {{R(1.14136170e-0), R(5.18118737e-1), R(1.85731975e-0)}},
    {{R(2.18708710e-1), R(1.93916682e-2), R(2.65755396e-0)}},
    {{R(1.68217399e-1), R(9.75705606e-1), R(7.24187871e-3)}},
    {{R(2.71719918e-1), R(4.19482500e-0), R(1.99325718e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(36, 0, // KR
    {{R(3.24386970e-1), R(1.31732163e-0), R(1.79912614e-0)}},
    {{R(6.31317973e-1), R(2.54706036e-1), R(3.23668394e-0)}},
    {{R(4.29961425e-3), R(1.00429433e-0), R(1.62188197e-1)}},
    {{R(1.98965610e-2), R(3.61094513e-0), R(2.45583672e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(37, 0, // RB
    {{R(2.90445351e-1), R(2.44201329e-0), R(7.69435449e-1)}},
    {{R(3.68420227e-2), R(1.16013332e-0), R(1.69591472e-1)}},
    {{R(1.58687000e-0), R(2.81617593e-3), R(1.28663830e-1)}},
    {{R(2.53082574e-0), R(1.88577417e-2), R(2.10753969e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(38, 0, // SR
    {{R(1.37373086e-2), R(1.97548672e-0), R(1.59261029e-0)}},
    {{R(1.87469061e-2), R(6.36079230e-0), R(2.21992482e-1)}},
    {{R(1.73263882e-1), R(4.66280378e-0), R(1.61265063e-3)}},
    {{R(2.01624958e-1), R(2.53027803e-1), R(1.53610568e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(39, 0, // Y
    {{R(6.75302747e-1), R(4.70286720e-1), R(2.63497677e-0)}},
    {{R(6.54331847e-2), R(1.06108709e-2), R(2.06643540e-0)}},
    {{R(1.09621746e-1), R(9.60348773e-1), R(5.28921555e-3)}},
    {{R(1.93131925e-1), R(1.63310938e-0), R(1.66083821e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(40, 0, // ZR
    {{R(2.64365505e-0), R(5.54225147e-1), R(7.61376625e-1)}},
    {{R(2.20202699e-0), R(1.78260107e-2), R(7.67218745e-2)}},
    {{R(6.02946891e-3), R(9.91630530e-2), R(9.56782020e-1)}},
    {{R(1.55143296e-2), R(1.76175995e-1), R(1.54330682e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(41, 0, // NB
    {{R(6.59532875e-1), R(1.84545854e-0), R(1.25584405e-0)}},
    {{R(8.66145490e-2), R(5.94774398e-0), R(6.40851475e-1)}},
    {{R(1.22253422e-1), R(7.06638328e-1), R(2.62381591e-3)}},
    {{R(1.66646050e-1), R(1.62853268e-0), R(8.26257859e-3)}})));
  addElement(ConstElementPtr(new ElementKirkland(42, 0, // MO
    {{R(6.10160120e-1), R(1.26544000e-0), R(1.97428762e-0)}},
    {{R(9.11628054e-2), R(5.06776025e-1), R(5.89590381e-0)}},
    {{R(6.48028962e-1), R(2.60380817e-3), R(1.13887493e-1)}},
    {{R(1.46634108e-0), R(7.84336311e-3), R(1.55114340e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(43, 0, // TC
    {{R(8.55189183e-1), R(1.66219641e-0), R(1.45575475e-0)}},
    {{R(1.02962151e-1), R(7.64907000e-0), R(1.01639987e-0)}},
    {{R(1.05445664e-1), R(7.71657112e-1), R(2.20992635e-3)}},
    {{R(1.42303338e-1), R(1.34659349e-0), R(7.90358976e-3)}})));
  addElement(ConstElementPtr(new ElementKirkland(44, 0, // RU
    {{R(4.70847093e-1), R(1.58180781e-0), R(2.02419818e-0)}},
    {{R(9.33029874e-2), R(4.52831347e-1), R(7.11489023e-0)}},
    {{R(1.97036257e-3), R(6.26912639e-1), R(1.02641320e-1)}},
    {{R(7.56181595e-3), R(1.25399858e-0), R(1.33786087e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(45, 0, // RH
    {{R(4.20051553e-1), R(1.76266507e-0), R(2.02735641e-0)}},
    {{R(9.38882628e-2), R(4.64441687e-1), R(8.19346046e-0)}},
    {{R(1.45487176e-3), R(6.22809600e-1), R(9.91529915e-2)}},
    {{R(7.82704517e-3), R(1.17194153e-0), R(1.24532839e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(46, 0, // PD
    {{R(2.10475155e-0), R(2.03884487e-0), R(1.82067264e-1)}},
    {{R(8.68606470e-0), R(3.78924449e-1), R(1.42921634e-1)}},
    {{R(9.52040948e-2), R(5.91445248e-1), R(1.13328676e-3)}},
    {{R(1.17125900e-1), R(1.07843808e-0), R(7.80252092e-3)}})));
  addElement(ConstElementPtr(new ElementKirkland(47, 0, // AG
    {{R(2.07981390e-0), R(4.43170726e-1), R(1.96515215e-0)}},
    {{R(9.92540297e-0), R(1.04920104e-1), R(6.40103839e-1)}},
    {{R(5.96130591e-1), R(4.78016333e-1), R(9.46458470e-2)}},
    {{R(8.89594790e-1), R(1.98509407e-0), R(1.12744464e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(48, 0, // CD
    {{R(1.63657549e-0), R(2.17927989e-0), R(7.71300690e-1)}},
    {{R(1.24540381e-1), R(1.45134660e-0), R(1.26695757e-1)}},
    {{R(6.64193880e-1), R(7.64563285e-1), R(8.61126689e-2)}},
    {{R(7.77659202e-1), R(1.66075210e-0), R(1.05728357e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(49, 0, // IN
    {{R(2.24820632e-0), R(1.64706864e-0), R(7.88679265e-1)}},
    {{R(1.51913507e-0), R(1.30113424e-1), R(1.06128184e-1)}},
    {{R(8.12579069e-2), R(6.68280346e-1), R(6.38467475e-1)}},
    {{R(9.94045620e-2), R(1.49742063e-0), R(7.18422635e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(50, 0, // SN
    {{R(2.16644620e-0), R(6.88691021e-1), R(1.92431751e-0)}},
    {{R(1.13174909e-1), R(1.10131285e-1), R(6.74464853e-1)}},
    {{R(5.65359888e-1), R(9.18683861e-1), R(7.80542213e-2)}},
    {{R(7.33564610e-1), R(1.02310312e-1), R(9.31104308e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(51, 0, // SB
    {{R(1.73662114e-0), R(9.99871380e-1), R(2.13972409e-0)}},
    {{R(8.84334719e-1), R(1.38462121e-1), R(1.19666432e-1)}},
    {{R(5.60566526e-1), R(9.93772747e-1), R(7.37374982e-2)}},
    {{R(6.72672880e-1), R(8.72330411e-0), R(8.78577715e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(52, 0, // TE
    {{R(2.09383882e-0), R(1.56940519e-0), R(1.30941993e-0)}},
    {{R(1.26856869e-1), R(1.21236537e-0), R(1.66633292e-1)}},
    {{R(6.98067804e-2), R(1.04969537e-0), R(5.55594354e-1)}},
    {{R(8.30817576e-2), R(7.43147857e-0), R(6.17487676e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(53, 0, // I
    {{R(1.60186925e-0), R(1.98510264e-0), R(1.48226200e-0)}},
    {{R(1.95031538e-1), R(1.36976183e-1), R(1.80304795e-0)}},
    {{R(5.53807199e-1), R(1.11728722e-0), R(6.60720847e-2)}},
    {{R(5.67912340e-1), R(6.40879878e-0), R(7.86615429e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(54, 0, // XE
    {{R(1.60015487e-0), R(1.71644581e-0), R(1.84968351e-0)}},
    {{R(2.92913354e-0), R(1.55882990e-1), R(2.22525983e-1)}},
    {{R(6.23813648e-2), R(1.21387555e-0), R(5.54051946e-1)}},
    {{R(7.45581223e-2), R(5.56013271e-0), R(5.21994521e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(55, 0, // CS
    {{R(2.95236854e-0), R(4.28105721e-1), R(1.89599233e-0)}},
    {{R(6.01461952e-0), R(4.64151246e-1), R(1.80109756e-1)}},
    {{R(5.48012938e-2), R(4.70838600e-0), R(5.90356719e-1)}},
    {{R(7.12799633e-2), R(4.56702799e-1), R(4.70236310e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(56, 0, // BA
    {{R(3.19434243e-0), R(1.98289586e-0), R(1.55121052e-1)}},
    {{R(9.27352241e-0), R(2.28741632e-1), R(3.82000231e-2)}},
    {{R(6.73222354e-2), R(4.48474211e-0), R(5.42674414e-1)}},
    {{R(7.30961745e-2), R(2.95703565e-1), R(4.08647015e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(57, 0, // LA
    {{R(2.05036425e-0), R(1.42114311e-1), R(3.23538151e-0)}},
    {{R(2.20348417e-1), R(3.96438056e-2), R(9.56979169e-0)}},
    {{R(6.34683429e-2), R(3.97960586e-0), R(5.20116711e-1)}},
    {{R(6.92443091e-2), R(2.53178406e-1), R(3.83614098e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(58, 0, // CE
    {{R(3.22990759e-0), R(1.57618307e-1), R(2.13477838e-0)}},
    {{R(9.94660135e-0), R(4.15378676e-2), R(2.40480572e-1)}},
    {{R(5.01907609e-1), R(3.80889010e-0), R(5.96625028e-2)}},
    {{R(3.66252019e-1), R(2.43275968e-1), R(6.59653503e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(59, 0, // PR
    {{R(1.58189324e-1), R(3.18141995e-0), R(2.27622140e-0)}},
    {{R(3.91309056e-2), R(1.04139545e-1), R(2.81671757e-1)}},
    {{R(3.97705472e-0), R(5.58448277e-2), R(4.85207954e-1)}},
    {{R(2.61872978e-1), R(6.30921695e-2), R(3.54234369e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(60, 0, // ND
    {{R(1.81379417e-1), R(3.17616396e-0), R(2.35221519e-0)}},
    {{R(4.37324793e-2), R(1.07842572e-1), R(3.05571833e-1)}},
    {{R(3.83125763e-0), R(5.25889976e-2), R(4.70090742e-1)}},
    {{R(2.54745408e-1), R(6.02676073e-2), R(3.39017003e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(61, 0, // PM
    {{R(1.92986811e-1), R(2.43756023e-0), R(3.17248504e-0)}},
    {{R(4.37785970e-2), R(3.29336996e-1), R(1.11259996e-1)}},
    {{R(3.58105414e-0), R(4.56529394e-1), R(4.94812177e-2)}},
    {{R(2.46709586e-1), R(3.24990282e-1), R(5.76553100e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(62, 0, // SM
    {{R(2.12002595e-1), R(3.16891754e-0), R(2.51503494e-0)}},
    {{R(4.57703608e-2), R(1.14536599e-1), R(3.55561054e-1)}},
    {{R(4.44080845e-1), R(3.36742101e-0), R(4.65652543e-2)}},
    {{R(3.11953363e-1), R(2.40291435e-1), R(5.52266819e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(63, 0, // EU
    {{R(2.59355002e-0), R(3.16557522e-0), R(2.29402652e-1)}},
    {{R(3.82452612e-1), R(1.17675155e-1), R(4.76642249e-2)}},
    {{R(4.32257780e-1), R(3.17261920e-0), R(4.37958317e-2)}},
    {{R(2.99719833e-1), R(2.34462738e-1), R(5.29440680e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(64, 0, // GD
    {{R(3.19144939e-0), R(2.55766431e-0), R(3.32681934e-1)}},
    {{R(1.20224655e-1), R(4.08338876e-1), R(5.85819814e-2)}},
    {{R(4.14243130e-2), R(2.61036728e-0), R(4.20526863e-1)}},
    {{R(5.06771477e-2), R(1.99344244e-1), R(2.85686240e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(65, 0, // TB
    {{R(2.59407462e-1), R(3.16177855e-0), R(2.75095751e-0)}},
    {{R(5.04689354e-2), R(1.23140183e-1), R(4.38337626e-1)}},
    {{R(2.79247686e-0), R(3.85931001e-2), R(4.10881708e-1)}},
    {{R(2.23797309e-1), R(4.87920992e-2), R(2.77622892e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(66, 0, // DY
    {{R(3.16055396e-0), R(2.82751709e-0), R(2.75140255e-1)}},
    {{R(1.25470414e-1), R(4.67899094e-1), R(5.23226982e-2)}},
    {{R(4.00967160e-1), R(2.63110834e-0), R(3.61333817e-2)}},
    {{R(2.67614884e-1), R(2.19498166e-1), R(4.68871497e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(67, 0, // HO
    {{R(2.88642467e-1), R(2.90567296e-0), R(3.15960159e-0)}},
    {{R(5.40507687e-2), R(4.97581077e-1), R(1.27599505e-1)}},
    {{R(3.91280259e-1), R(2.48596038e-0), R(3.37664478e-2)}},
    {{R(2.58151831e-1), R(2.15400972e-1), R(4.50664323e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(68, 0, // ER
    {{R(3.15573213e-0), R(3.11519560e-1), R(2.97722406e-0)}},
    {{R(1.29729009e-1), R(5.81399387e-2), R(5.31213394e-1)}},
    {{R(3.81563854e-1), R(2.40247532e-0), R(3.15224214e-2)}},
    {{R(2.49195776e-1), R(2.13627616e-1), R(4.33253257e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(69, 0, // TM
    {{R(3.15591970e-0), R(3.22544710e-1), R(3.05569053e-0)}},
    {{R(1.31232407e-1), R(5.97223323e-2), R(5.61876773e-1)}},
    {{R(2.92845100e-2), R(3.72487205e-1), R(2.27833695e-0)}},
    {{R(4.16534255e-2), R(2.40821967e-1), R(2.10034185e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(70, 0, // YB
    {{R(3.10794704e-0), R(3.14091221e-0), R(3.75660454e-1)}},
    {{R(6.06347847e-1), R(1.33705269e-1), R(7.29814740e-2)}},
    {{R(3.61901097e-1), R(2.45409082e-0), R(2.72383990e-2)}},
    {{R(2.32652051e-1), R(2.12695209e-1), R(3.99969597e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(71, 0, // LU
    {{R(3.11446863e-0), R(5.39634353e-1), R(3.06460915e-0)}},
    {{R(1.38968881e-1), R(8.91708508e-2), R(6.79919563e-1)}},
    {{R(2.58563745e-2), R(2.13983556e-0), R(3.47788231e-1)}},
    {{R(3.82808522e-2), R(1.80078788e-1), R(2.22706591e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(72, 0, // HF
    {{R(3.01166899e-0), R(3.16284788e-0), R(6.33421771e-1)}},
    {{R(7.10401889e-1), R(1.38262192e-1), R(9.48486572e-2)}},
    {{R(3.41417198e-1), R(1.53566013e-0), R(2.40723773e-2)}},
    {{R(2.14129678e-1), R(1.55298698e-1), R(3.67833690e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(73, 0, // TA
    {{R(3.20236821e-0), R(8.30098413e-1), R(2.86552297e-0)}},
    {{R(1.38446369e-1), R(1.18381581e-1), R(7.66369118e-1)}},
    {{R(2.24813887e-2), R(1.40165263e-0), R(3.33740596e-1)}},
    {{R(3.52934622e-2), R(1.46148877e-1), R(2.05704486e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(74, 0, // W
    {{R(9.24906855e-1), R(2.75554557e-0), R(3.30440060e-0)}},
    {{R(1.28663377e-1), R(7.65826479e-1), R(1.34471170e-1)}},
    {{R(3.29973862e-1), R(1.09916444e-0), R(2.06498883e-2)}},
    {{R(1.98218895e-1), R(1.35087534e-1), R(3.38918459e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(75, 0, // RE
    {{R(1.96952105e-0), R(1.21726619e-0), R(4.10391685e-0)}},
    {{R(4.98830620e-1), R(1.33243809e-1), R(1.84396916e-0)}},
    {{R(2.90791978e-2), R(2.30696669e-1), R(6.08840299e-1)}},
    {{R(2.84192813e-2), R(1.90968784e-1), R(1.37090356e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(76, 0, // OS
    {{R(2.06385867e-0), R(1.29603406e-0), R(3.96920673e-0)}},
    {{R(4.05671697e-1), R(1.46559047e-1), R(1.82561596e-0)}},
    {{R(2.69835487e-2), R(2.31083999e-1), R(6.30466774e-1)}},
    {{R(2.84172045e-2), R(1.79765184e-1), R(1.38911543e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(77, 0, // IR
    {{R(2.21522726e-0), R(1.37573155e-0), R(3.78244405e-0)}},
    {{R(3.24464090e-1), R(1.60920048e-1), R(1.78756553e-0)}},
    {{R(2.44643240e-2), R(2.36932016e-1), R(6.48471412e-1)}},
    {{R(2.82909938e-2), R(1.70692368e-1), R(1.37928390e-0)}})));
  addElement(ConstElementPtr(new ElementKirkland(78, 0, // PT
    {{R(9.84697940e-1), R(2.73987079e-0), R(3.61696715e-0)}},
    {{R(1.60910839e-1), R(7.18971667e-1), R(1.29281016e-1)}},
    {{R(3.02885602e-1), R(2.78370726e-1), R(1.52124129e-2)}},
    {{R(1.70134854e-1), R(1.49862703e-0), R(2.83510822e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(79, 0, // AU
    {{R(9.61263398e-1), R(3.69581030e-0), R(2.77567491e-0)}},
    {{R(1.70932277e-1), R(1.29335319e-1), R(6.89997070e-1)}},
    {{R(2.95414176e-1), R(3.11475743e-1), R(1.43237267e-2)}},
    {{R(1.63525510e-1), R(1.39200901e-0), R(2.71265337e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(80, 0, // HG
    {{R(1.29200491e-0), R(2.75161478e-0), R(3.49387949e-0)}},
    {{R(1.83432865e-1), R(9.42368371e-1), R(1.46235654e-1)}},
    {{R(2.77304636e-1), R(4.30232810e-1), R(1.48294351e-2)}},
    {{R(1.55110144e-1), R(1.28871670e-0), R(2.61903834e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(81, 0, // TL
    {{R(3.75964730e-0), R(3.21195904e-0), R(6.47767825e-1)}},
    {{R(1.35041513e-1), R(6.66330993e-1), R(9.22518234e-2)}},
    {{R(2.76123274e-1), R(3.18838810e-1), R(1.31668419e-2)}},
    {{R(1.50312897e-1), R(1.12565588e-0), R(2.48879842e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(82, 0, // PB
    {{R(1.00795975e-0), R(3.09796153e-0), R(3.61296864e-0)}},
    {{R(1.17268427e-1), R(8.80453235e-1), R(1.47325812e-1)}},
    {{R(2.62401476e-1), R(4.05621995e-1), R(1.31812509e-2)}},
    {{R(1.43491014e-1), R(1.04103506e-0), R(2.39575415e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(83, 0, // BI
    {{R(1.59826875e-0), R(4.38233925e-0), R(2.06074719e-0)}},
    {{R(1.56897471e-1), R(2.47094692e-0), R(5.72438972e-1)}},
    {{R(1.94426023e-1), R(8.22704978e-1), R(2.33226953e-2)}},
    {{R(1.32979109e-1), R(9.56532528e-1), R(2.23038435e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(84, 0, // PO
    {{R(1.71463223e-0), R(2.14115960e-0), R(4.37512413e-0)}},
    {{R(9.79262841e-1), R(2.10193717e-1), R(3.66948812e-0)}},
    {{R(2.16216680e-2), R(1.97843837e-1), R(6.52047920e-1)}},
    {{R(1.98456144e-2), R(1.33758807e-1), R(7.80432104e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(85, 0, // AT
    {{R(1.48047794e-0), R(2.09174630e-0), R(4.75246033e-0)}},
    {{R(1.25943919e-2), R(1.83803008e-1), R(4.19890596e-0)}},
    {{R(1.85643958e-2), R(2.05859375e-1), R(7.13540948e-1)}},
    {{R(1.81383503e-2), R(1.33035404e-1), R(7.03031938e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(86, 0, // RN
    {{R(6.30022295e-1), R(3.80962881e-0), R(3.89756067e-0)}},
    {{R(1.40909762e-1), R(3.08515540e-1), R(6.51559763e-1)}},
    {{R(2.40755100e-1), R(2.62868577e-0), R(3.14285931e-2)}},
    {{R(1.08899672e-1), R(6.42383261e-0), R(2.42346699e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(87, 0, // FR
    {{R(5.23288135e-0), R(2.48604205e-0), R(3.23431354e-1)}},
    {{R(8.60599536e-0), R(3.04543982e-1), R(3.87759096e-2)}},
    {{R(2.55403596e-1), R(5.53607228e-1), R(5.75278889e-3)}},
    {{R(1.28717724e-1), R(5.36977452e-1), R(1.29417790e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(88, 0, // RA
    {{R(1.44192685e-0), R(3.55291725e-0), R(3.91259586e-0)}},
    {{R(1.18740873e-1), R(1.01739750e-0), R(6.31814783e-1)}},
    {{R(2.16173519e-1), R(3.94191605e-0), R(4.60422605e-2)}},
    {{R(9.55806441e-2), R(3.50602732e-1), R(2.20850385e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(89, 0, // AC
    {{R(1.45864127e-0), R(4.18945405e-0), R(3.65866182e-0)}},
    {{R(1.07760494e-1), R(8.89090649e-1), R(1.05088931e-0)}},
    {{R(2.08479229e-1), R(3.16528117e-0), R(5.23892556e-2)}},
    {{R(9.09335557e-2), R(3.13297788e-1), R(2.08807697e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(90, 0, // TH
    {{R(1.19014064e-0), R(2.55380607e-0), R(4.68110181e-0)}},
    {{R(7.73468729e-2), R(6.59693681e-1), R(1.28013896e-1)}},
    {{R(2.26121303e-1), R(3.58250545e-1), R(7.82263950e-3)}},
    {{R(1.08632194e-1), R(4.56765664e-1), R(1.62623474e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(91, 0, // PA
    {{R(4.68537504e-0), R(2.98413708e-0), R(8.91988061e-1)}},
    {{R(1.44503632e-1), R(5.56438592e-1), R(6.69512914e-2)}},
    {{R(2.24825384e-1), R(3.04444846e-1), R(9.48162708e-3)}},
    {{R(1.03235396e-1), R(4.27255647e-1), R(1.77730611e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(92, 0, // U
    {{R(4.63343606e-0), R(3.18157056e-0), R(8.76455075e-1)}},
    {{R(1.63377267e-1), R(5.69517868e-1), R(6.88860012e-2)}},
    {{R(2.21685477e-1), R(2.72917100e-1), R(1.11737298e-2)}},
    {{R(9.84254550e-2), R(4.09470917e-1), R(1.86215410e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(93, 0, // NP
    {{R(4.56773888e-0), R(3.40325179e-0), R(8.61841923e-1)}},
    {{R(1.90992795e-1), R(5.90099634e-1), R(7.03204851e-2)}},
    {{R(2.19728870e-1), R(2.38176903e-1), R(1.38306499e-2)}},
    {{R(9.36334280e-2), R(3.93554882e-1), R(1.94437286e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(94, 0, // PU
    {{R(5.45671123e-0), R(1.11687906e-1), R(3.30260343e-0)}},
    {{R(1.01892720e-1), R(3.98131313e-2), R(3.14622212e-1)}},
    {{R(1.84568319e-1), R(4.93644263e-1), R(3.57484743e-0)}},
    {{R(1.04220860e-1), R(4.63080540e-1), R(2.19369542e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(95, 0, // AM
    {{R(5.38321999e-0), R(1.23343236e-1), R(3.46469090e-0)}},
    {{R(1.07289857e-1), R(4.15137806e-2), R(3.39326208e-1)}},
    {{R(1.75437132e-1), R(3.39800073e-0), R(4.69459519e-1)}},
    {{R(9.98932346e-2), R(2.11601535e-1), R(4.51996970e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(96, 0, // CM
    {{R(5.38402377e-0), R(3.49861264e-0), R(1.88039547e-1)}},
    {{R(1.11211419e-1), R(3.56750210e-1), R(5.39853583e-2)}},
    {{R(1.69143137e-1), R(3.19595016e-0), R(4.64393059e-1)}},
    {{R(9.60082633e-2), R(1.80694389e-1), R(4.36318197e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(97, 0, // BK
    {{R(3.66090688e-0), R(2.03054678e-1), R(5.30697515e-0)}},
    {{R(3.84420906e-1), R(5.48547131e-2), R(1.17150262e-1)}},
    {{R(1.60934046e-1), R(3.04808401e-0), R(4.43610295e-1)}},
    {{R(9.21020329e-2), R(1.73525367e-1), R(4.27132359e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(98, 0, // CF
    {{R(3.94150390e-0), R(5.16915345e-0), R(1.61941074e-1)}},
    {{R(4.18246722e-1), R(1.25201788e-1), R(4.81540117e-2)}},
    {{R(4.15299561e-1), R(2.91761325e-0), R(1.51474927e-1)}},
    {{R(4.24913856e-1), R(1.90899693e-1), R(8.81568925e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(99, 0, // ES
    {{R(4.09780623e-0), R(5.10079393e-0), R(1.74617289e-1)}},
    {{R(4.46021145e-1), R(1.31768613e-1), R(5.02742829e-2)}},
    {{R(2.76774658e-0), R(1.44496639e-1), R(4.02772109e-1)}},
    {{R(1.84815393e-1), R(8.46232592e-2), R(4.17640100e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(100, 0, // FM
    {{R(4.24934820e-0), R(5.03556594e-0), R(1.88920613e-1)}},
    {{R(4.75263933e-1), R(1.38570834e-1), R(5.26975158e-2)}},
    {{R(3.94356058e-1), R(2.61213100e-0), R(1.38001927e-1)}},
    {{R(4.11193751e-1), R(1.78537905e-1), R(8.12774434e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(101, 0, // MD
    {{R(2.00942931e-1), R(4.40119869e-0), R(4.97250102e-0)}},
    {{R(5.48366518e-2), R(5.04248434e-1), R(1.45721366e-1)}},
    {{R(2.47530599e-0), R(3.86883197e-1), R(1.31936095e-1)}},
    {{R(1.72978308e-1), R(4.05043898e-1), R(7.80821071e-2)}})));
  addElement(ConstElementPtr(new ElementKirkland(102, 0, // NO
    {{R(2.16052899e-1), R(4.91106799e-0), R(4.54862870e-0)}},
    {{R(5.83584058e-2), R(1.53264212e-1), R(5.34434760e-1)}},
    {{R(2.36114249e-0), R(1.26277292e-1), R(3.81364501e-1)}},
    {{R(1.68164803e-1), R(7.50304633e-2), R(3.99305852e-1)}})));
  addElement(ConstElementPtr(new ElementKirkland(103, 0, // LR
    {{R(4.86738014e-0), R(3.19974401e-1), R(4.58872425e-0)}},
    {{R(1.60320520e-1), R(6.70871138e-2), R(5.77039373e-1)}},
    {{R(1.21482448e-1), R(2.31639872e-0), R(3.79258137e-1)}},
    {{R(7.22275899e-2), R(1.41279737e-1), R(3.89973484e-1)}})));
}

class AURORA_API RegisterKirkland
{
public:
  RegisterKirkland()
  {
    FormFactorParametrization::registerClass<FormFactorParametrizationKirkland>
      ("Kirkland");
  }
} registerKirkland;

} // namespace Aurora
