//--- Aurora/AuroraLib/Energy.hpp ----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_ENERGY_HPP
#define AURORA_AURORA_LIB_ENERGY_HPP

#include "AuroraLib/Physics.hpp"

namespace Aurora
{
class AURORA_API Energy
{
public:
  Energy() = default;

  Energy(Real energy)
  {
    setEnergy(energy);
  }

  /// Returns the energy of the electron in eV
  Real energy() const
  {
    return mEnergy;
  }

  operator Real() const
  {
    return mEnergy;
  }

  /// Sets the kinetic energy of the electron in eV
  void setEnergy(Real energy)
  {
    using namespace Phys;

    assert(energy >= R(0.0));
    mEnergy = energy;

    #if (1 == AURORA_NO_RELATIVISTICS)
      mWavenumber = Math::twoPi * hBar / sqrt(2 * mE * kineticEnergy);
    #else
      using namespace Phys;
      constexpr Real const1 = R(0.5) * mE / (Math::piSqr * hBar * hBar);
      constexpr Real const2 = R(0.25) / (Math::piSqr * hBar * hBar);
      mWavelength = R(1.0) / std::sqrt(const1 * mEnergy +
                                       const2 * mEnergy * mEnergy);
    #endif
  }

  /// Return the wavelength in angstrom. This function respects relativistic
  /// effects.
  /// @details
  ///  It uses the formula
  ///    @f[
  ///      \lambda = \frac{2\,\pi\,\hbar c} {\sqrt{2\,m_e\,c^2\,T+T^2}}
  ///    @f]
  Real wavelength() const
  {
    return mWavelength;
  }

  /// Returns the relativistic wavenumber @f$k@f$ of the electron in Angstrom.
  /// @details
  ///  It uses the formula
  ///  @f[
  ///    k = \frac{2\,\pi}{\lambda} = \frac{\sqrt{2\,m_e\,c^2\,T+T^2}}{\hbar\,c}
  ///  @f]
  ///  where @f$T@f$ is the kinetic energy of the electron.
  Real wavenumber() const
  {
    return Math::twoPi / mWavelength;
  }

  /// Returns the Lorentz factor or relativistic factor.
  /// @details
  ///  The function uses the formula
  ///    @f[
  ///      \gamma = 1 + \frac{T}{m_e c^2}
  ///    @f]
  Real lorentzFactor() const
  {
    #if (1 == AURORA_NO_RELATIVISTICS)
      return R(1.0);
    #else
      return R(1.0) + mEnergy / Phys::mE;
    #endif
  }

  /// Returns the speed of the electron in fractions of the speed of light
  Real speed() const
  {
    return lorentzFactorToSpeed(lorentzFactor());
  }

  bool operator==(const Energy& other) const
  {
    return (mEnergy == other.mEnergy);
  }

  bool operator!=(const Energy& other) const
  {
    return !operator==(other);
  }

private:
  // the energy of the electron in eV
  Real mEnergy = 0;
  Real mWavelength = 0;
};

inline std::istream& operator>>(std::istream& inStream, Energy& energy)
{
  Real temp;
  inStream >> temp;
  energy.setEnergy(temp);
  return inStream;
}

} // namespace Aurora

#endif
