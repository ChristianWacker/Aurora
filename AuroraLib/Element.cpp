//--- Aurora/AuroraLib/Element.cpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "Element.hpp"

#include "Exception.hpp"
#include "Math.hpp"
#include "Physics.hpp"
#include "Quadrature.hpp"

#include <cmath>
#include <limits>

namespace Aurora
{

/* static */ const size_t Element::mNumElements;

const std::string Element::mSymbols[Element::mNumElements] =
  {"H" , "HE", "LI", "BE", "B" , "C" , "N" , "O" , "F" , "NE",
   "NA", "MG", "AL", "SI", "P" , "S" , "CL", "AR", "K" , "CA",
   "SC", "TI", "V" , "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
   "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y" , "ZR",
   "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
   "SB", "TE", "I" , "XE", "CS", "BA", "LA", "CE", "PR", "ND",
   "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
   "LU", "HF", "TA", "W" , "RE", "OS", "IR", "PT", "AU", "HG",
   "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
   "PA", "U" , "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
   "MD", "NO", "LR"};

Element::Element(size_t atomicNumber, int charge) :
  mAtomicNumber(atomicNumber), mCharge(charge)
{
  assert((atomicNumber >= 1) && (atomicNumber <= mNumElements));
}

Real Element::phaseShift(Real rho, Real z, Real bFactor) const
{
  using namespace Phys;
  const Real rhoSqr = rho * rho;

  auto f = [rhoSqr, bFactor, this] (Real zPrime)
  {
    return potential(rhoSqr + zPrime * zPrime, bFactor);
  };

  // integrate using the double exponential method, as the integrand might
  // be singular for (rho = 0 = z)
  return -mE / (hBar * hBar) * integrateDE(f, Real(0.0), z);
}

Real Element::projectedPotential(Real, Real) const
{
  return std::numeric_limits<Real>::quiet_NaN();
}

Real Element::xRayFormFactor(Real) const
{
  return std::numeric_limits<Real>::quiet_NaN();
}

Real Element::hybridPotential(Real, Real, Real) const
{
  return std::numeric_limits<Real>::quiet_NaN();
}

Real Element::hybridPhaseShift(Real, Real, Real) const
{
 return std::numeric_limits<Real>::quiet_NaN();
}

Real Element::radiusSqr() const
{
  return std::numeric_limits<Real>::quiet_NaN();
}

std::string Element::symbol(size_t atomicNumber)
{
  assert((atomicNumber >= 1) && (atomicNumber <= mNumElements));
  return mSymbols[atomicNumber - 1];
}

size_t Element::atomicNumber(const std::string& symbol)
{
  for (size_t i = 1; i <= mNumElements; ++i)
    if (mSymbols[i - 1] == symbol)
      return i;

  return 0;
}

const Element& FormFactorParametrization::element(size_t atomicNumber,
                                                  int charge) const
{
  assert((atomicNumber >= 1) && (atomicNumber <= Element::mNumElements));
  auto iterator = mElementMap.find(std::make_pair(atomicNumber, charge));

  if (mElementMap.end() == iterator)
  {
    AURORA_THROW(EInvalidParameter,
                 "Element \"" + std::to_string(atomicNumber) + "\" not found.");
  }

  return *(iterator->second);
}

const Element& FormFactorParametrization::element(const std::string& symbol,
                                                  int charge) const
{
  return element(Element::atomicNumber(symbol), charge);
}

const Element& FormFactorParametrization::element(
  const std::string& extendedSymbol) const
{
  const size_t length = extendedSymbol.length();
  if (length > 2 ||
      '+' == extendedSymbol[length - 1] || '-' == extendedSymbol[length - 1])
  {
    if (isdigit(extendedSymbol[length - 2]))
    {
      int32_t charge = extendedSymbol[length - 2] - '0';
      if (extendedSymbol[length - 1] == '-')
        charge = -charge;

      const std::string symbol = extendedSymbol.substr(0, length - 2);
      return element(symbol, charge);
    }
    else
    {
      AURORA_THROW(EInvalidParameter, "Missing Digit.");
    }
  }

  return element(Element::atomicNumber(extendedSymbol), 0);
}

void FormFactorParametrization::addElement(ConstElementPtr element)
{
  const size_t atomicNumber = element->atomicNumber();
  const int32_t charge = element->charge();
  auto key = std::make_pair(atomicNumber, charge);

  if (!mElementMap.insert(std::make_pair(key, std::move(element))).second)
  {
    AURORA_THROW(EInvalidParameter, "Element \"" + element->symbol() +
                                    "\" has already been added.");
  }
}

FormFactorParametrization::Creators& FormFactorParametrization::creators()
{
  static Creators creators;
  return creators;
}

} // namespace Aurora
