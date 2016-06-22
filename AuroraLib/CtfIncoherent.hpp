//--- Aurora/AuroraLib/CtfIncoherent.hpp -------------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_CTF_INCOHERENT_HPP
#define AURORA_AURORA_LIB_CTF_INCOHERENT_HPP

#include "AuroraLib/Ctf.hpp"

namespace Aurora
{

class AURORA_API CtfIncoherent
{
public:
  CtfIncoherent(Real energy);

  CtfCoherent& coherent()
  {
    return mCoherent;
  }

  /// This function returns the incoherent CTF, i.e. the CTF including the
  /// damping envelopes.
  /// @param info
  ///  describes the properties of the target buffer.
  /// @param factor
  ///  Scaling factor for the whole CTF. This can be used to correct the FFT
  ///  scaling.
  CBuffer2D buffer(const Buffer2DInfoBase& info, Real factor = 1) const;

  /// @{
  /// Returns the standard deviation of the focus spread
  /// @f$\sigma(D)@f$
  /// @f{eqnarray*}{
  ///  \sigma(D)^2 &=& C_C^2 \left[\left(\frac{\Delta U}{U_0}\right)^2 +
  ///  4 \left(\frac{\Delta I}{I_0}\right)^2 + \left(\frac{\sigma(E)}{E_0}
  ///  \right)^2\right]
  /// @f}
  Real focusSpread() const
  {
    return mCC * std::sqrt(powerOf<2>(mHighTensionRipple) +
                           4 * powerOf<2>(mCurrentInstability) +
                           powerOf<2>(mEnergySpread / energy()));
  }

  /// Returns the temporal envelope function
  /// @f{eqnarray*}{
  ///  E_T(\vartheta) &=& \exp\left(-\frac{\sigma(D)^2\,k^2\,\vartheta^4}{8}
  ///    \right)
  /// @f}
  Real temporalEnvelope(Real theta) const;

  /// Returns the spatial envelope function
  /// @f{eqnarray*}{
  ///  E_S(\vartheta) &=& \exp\left(-\frac{\sigma(\beta)^2\,k^2\left(C_1\,
  ///    \vartheta + C_3\,\vartheta^3\right)^2}{2}\right)
  /// @f}
  Real spatialEnvelope(Real theta) const;
  /// @}

  /// Returns the coefficient of the chromatic aberration @f$C_C@f$ in Angstrom
  Real CC() const
  {
    return mCC;
  }

  /// Sets the coefficient of the chromatic aberration @f$C_C@f$ in Angstrom
  void setCC(Real cc)
  {
    mCC = cc;
  }

  /// Returns the high tension ripple @f$\frac{\Delta U}{U_0}@f$
  Real highTensionRipple() const
  {
    return mHighTensionRipple;
  }

  /// Sets the high tension ripple @f$\frac{\Delta U}{U_0}@f$
  void setHighTensionRipple(Real highTensionRipple)
  {
    mHighTensionRipple = highTensionRipple;
  }

  /// Returns the lens current instabilities @f$\frac{\Delta I}{I_0}@f$
  Real currentInstabilities() const
  {
    return mCurrentInstability;
  }

  /// Sets the lens current instabilities @f$\frac{\Delta I}{I_0}@f$
  void setCurrentInstability(Real currentInstability)
  {
    mCurrentInstability = currentInstability;
  }

  /// Return the energy spread @f$\sigma(E)@f$ of the electron gun in eV
  Real energySpread() const
  {
    return mEnergySpread;
  }

  /// Set the energy spread @f$\sigma(E)@f$ of the electron gun in eV
  void setEnergySpread(Real energySpread)
  {
    mEnergySpread = energySpread;
  }

  /// Returns the illumination semi-angle
  /// @f$\alpha_i=\sigma(\beta)\sqrt{2 \ln 2}@f$
  Real illuminationSemiAngle() const
  {
    return mSigmaBeta * std::sqrt(2 * std::log(R(2.0)));
  }

  /// Sets the illumination semi-angle
  /// @f$\alpha_i=\sigma(\beta)\sqrt{2 \ln 2}@f$
  void setIlluminationSemiAngle(Real semiAngle)
  {
    mSigmaBeta = semiAngle / std::sqrt(2 * std::log(R(2.0)));
  }

  /// Returns the spread of the incident angle
  Real sigmaBeta() const
  {
    return mSigmaBeta;
  }

  Real energy() const
  {
    return mCoherent.energy();
  }

  /// Set the energy in eV.
  void setEnergy(Real energy)
  {
    mCoherent.setEnergy(energy);
  }

  Real wavenumber() const
  {
    return mCoherent.wavenumber();
  }

private:
  CtfCoherent mCoherent;

  // the coefficient of the chromatic aberration
  Real mCC;
  Real mHighTensionRipple;
  Real mCurrentInstability;
  // the energy spread of the electron gun in eV
  Real mEnergySpread;

  // spread of the incident angle in rad
  Real mSigmaBeta;
};

} // namespace Aurora

#endif
