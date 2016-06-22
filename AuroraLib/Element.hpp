//--- Aurora/AuroraLib/Element.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_ELEMENT_HPP
#define AURORA_AURORA_LIB_ELEMENT_HPP

#include "AuroraLib/Factory.hpp"

#include <array>

namespace Aurora
{

/// Class to describe the physical properties of a chemical element, most
/// importantly it provides a way to calculate the phase shift.
class AURORA_API Element
{
public:
  virtual ~Element() { }

  /// Prevent copying
  Element(const Element&) = delete;
  Element& operator=(const Element&) = delete;

  /// Returns the phase shift caused by the element in 1 / Angstrom
  /// @f[
  ///   \Phi(\rho,\,z)=-\frac{m}{\hbar^2} \intop_0^zV(\rho,\,z')\text{d}z',
  /// @f].
  /// @remarks
  ///  The calculated phase shift includes a Debye-Waller factor that is
  ///  specified by @p bFactor
  ///  @f[
  ///   \text{DW} = \exp(-B\,K^2)
  ///  @f].
  virtual Real phaseShift(Real rho, Real z, Real bFactor = Real()) const;

  /// Returns the potential projected along one axis
  ///  @f[
  ///    \intop_{-\infty}^{\infty}V(\rho,\,z') \text{d}z'
  ///  @f]
  /// @param rhoSqr is the squared distance in Angstrom^2 of the projection line
  ///  and the atom.
  /// @param bFactor specifies the Debye-Waller factor see
  ///  electronicFormFactor()
  virtual Real projectedPotential(Real rhoSqr, Real bFactor = 0) const;

  /// Return the electronic form factor
  ///  @f[
  ///  f_{e}\left(\mathbf{K}\right)=-\frac{m}{2\,\pi\,\hbar^{2}}\int \
  ///    \exp\left(\text{i}\,\mathbf{K}\cdot\mathbf{r}'\right)\
  ///    V\left(\mathbf{r}'\right)\text{d}^{3}r',
  ///  @f]
  /// The form factor is multiplied by a Debye-Waller (DW) factor that is
  /// specified by @p bFactor
  ///  @f[
  ///   \text{DW} = \exp(-B\,K^2)
  ///  @f]
  /// where @f$B@f$ is connected to the mean square displacement
  /// @f$\left<u^2\right>@f$ by
  ///  @f[
  ///   B = 8\,\pi^2\,\left<u^2\right>
  ///  @f]
  /// @param kSqr
  ///  must be in in Anstrom^-2
  /// @param bFactor
  ///  must be in in Anstrom^-2
  virtual Real electronicFormFactor(Real kSqr, Real bFactor = 0) const = 0;

  /// Returns the X-ray form factor
  ///  @f[
  ///  f_{\gamma}\left(\mathbf{K}\right)=-\int \
  ///    \exp\left(\text{i}\,\mathbf{K}\cdot\mathbf{r}'\right)\
  ///    \rho_e\left(\mathbf{r}'\right)\text{d}^{3}r',
  ///  @f]
  /// where @f$\rho_E(r)@f$ denotes the electron charge density.
  /// @retval
  ///  The return value has no physical dimension.
  /// @param kSqr
  ///  must be in in Anstrom^-2
  virtual Real xRayFormFactor(Real kSqr) const;

  /// Return the relativistically uncorrected potential energy.
  /// @retval
  ///  The energy is returned in eV.
  /// @param rSqr
  ///  must be in Angstrom^2.
  /// @param bFactor
  ///  must be in Angstrom^2.
  virtual Real potential(Real rSqr, Real bFactor = 0) const = 0;

  virtual Real hybridPotential(Real kappaSqr, Real z, Real bFactor) const;

  virtual Real hybridPhaseShift(Real kappaSqr, Real z, Real bFactor) const;

  /// Return the mean quadratic distance of the electron distribution
  /// @retval
  ///  The distance is returned in Angstrom^2.
  virtual Real radiusSqr() const;

  /// Return the symbol of the element.
  std::string symbol() const
  {
    return mSymbols[mAtomicNumber - 1];
  }

  /// Returns the atomic number of the element.
  size_t atomicNumber() const
  {
    return mAtomicNumber;
  }

  int charge() const
  {
    return mCharge;
  }

  /// Returns the the symbol of chemical element given the atomic number.
  static std::string symbol(size_t atomicNumber);
  /// Returns the atomic number of the element given symbol. If the symbol
  /// can't be found 0 is returned.
  static size_t atomicNumber(const std::string& symbol);

  // The number of chemical elements, for which the symbols are known.
  static const size_t mNumElements = 103;

protected:
  /// Creates a new element.
  /// @param atomicNumber
  ///  of the chemical element, e.g. Z = 14 for silicon.
  /// @param charge
  ///  of the chemical elements in units of the elementary charge. This is used
  ///  for the simulation of ions.
  Element(size_t atomicNumber, int charge);

  size_t mAtomicNumber;
  int mCharge;

private:
  // The symbols of the chemical elements.
  static const std::string mSymbols[mNumElements];
};
typedef std::unique_ptr<const Element> ConstElementPtr;

class AURORA_API FormFactorParametrization :
  public Factory<FormFactorParametrization>
{
public:
  const Element& element(size_t atomicNumber, int charge) const;
  const Element& element(const std::string& symbol, int charge) const;
  const Element& element(const std::string& extendedSymbol) const;

  static Creators& creators();

protected:
  void addElement(ConstElementPtr element);

private:
#if (AURORA_COMPILER == AURORA_COMPILER_MSVC)
  std::map<std::pair<size_t, int>,
           std::shared_ptr<const Element>> mElementMap;
#else
  std::map<std::pair<size_t, int>,
           std::unique_ptr<const Element>> mElementMap;
#endif
};

} // namespace Aurora

#endif
