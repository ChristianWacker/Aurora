//--- Aurora/AuroraLib/Ctf.hpp -------------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_CTF_HPP
#define AURORA_AURORA_LIB_CTF_HPP

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/Energy.hpp"
#include "AuroraLib/Math.hpp"

#include <array>

namespace Aurora
{

/// @brief
///  Class for describing the contrast transfer function (CTF) of a
///  transmission electron microscope. It models the effect of the imaging
///  system on the exit-wavefunction, i.e. the wavefunction of the electron
///  after the sample.
/// @details
///   This class includes all axial aberrations up to and including the
///   seventh order. We use the notation of [Haider2008] and [Erni2010]. The
///   phase shift @f$\chi@f$ is given by [Erni2010], p. 221, eq. 7.14:
///   @f{eqnarray*}{
///     \chi(\omega) &=& \text{Re}\left(A_0\,\omega^* + C_1\frac{\omega\,
///     \omega^*}{2} + A_1\frac{\omega^{*2}}{2} + B_2\omega^2\,\omega^* +
///     A_2\frac{\omega^{*3}}{3}\right.\\
///     &&+ S_3 \omega^3 \omega^* + C_3 \frac{(\omega\omega^*)^2}{4} +
///     A_3 \frac{\omega^{*4}}{4}\\
///     &&+ D_4 \omega^4 \omega^* + B_4 \omega^3 \omega^{*2} +
///     A_4 \frac{\omega^{*5}}{5}\\
///     &&+ R_5 \omega^5 \omega^* + S_5 \omega^4 \omega^{*2} +
///     C_5 \frac{(\omega \omega^*)^3}{6} + A_5 \frac{\omega^{6*}}{6} \\
///     &&+ F_6 \omega^6 \omega^* + D_6 \omega^5 \omega^{*2} +
///     B_6 \omega^4 \omega^{*3} + A_6 \frac{\omega^{*7}}{7}\\
///     &&+ \left. G_7 \omega^7 \omega^* + R_7 \omega^6 \omega^{*2} +
///     S_7 \omega^5 \omega^{*3} + C_7 \frac{(\omega \omega^*)^4}{8} + A_7
///     \frac{\omega^{* 8}}{8} \right)
///   @f}
///   where @f$\omega=\vartheta_x+\text{i}\,\vartheta_y@f$ describes a point
///   in the backfocal plane. @f$\vartheta_x@f$ and @f$\vartheta_x@f$ are the
///   scattering angles in the respective directions.
///   The numerical subscript of each aberrations indicates its order @f$n@f$:
///    - @f$A_n@f$ is the \f$n\f$-th order astigmatism which has an
///      @f$(n+1)@f$-fold symmetry.
///    - @f$C_n@f$ is the @f$n@f$-th order spherical aberration, which
///      is isotropic. @f$C_1@f$ is the defocus and @f$C_3@f$ is commonly
///      known as @f$C_s@f$.
///    - @f$B_n@f$ is the @f$n@f$-th order axial coma of symmetry 1.
///    - @f$S_n@f$ is the @f$n@f$-th order star aberration with twofold
///      symmetry.
///    - @f$D_n@f$ is the @f$n@f$-th order three-lobe aberration with
///      threefold symmetry.
///    - @f$R_n@f$ is the @f$n@f$-th order rosette aberration with fourfold
///      symmetry.
///    - @f$F_6@f$ is the sixth-order pentacle aberration with fivefold
///      symmetry.
///    - @f$G_7@f$ is the seventh-order chaplet aberration with sixfold
///      symmetry.
///
///   With the exception of the isotropic, spherical coefficients @f$C_n@f$,
///   all the aberrations coefficients are complex-valued.
class AURORA_API CtfCoherent : public Energy
{
public:
  using Energy::Energy;

  /// Returns the phase shift @f$k\,\chi@f$ for the given scattering angles
  /// @param defocus
  ///  Additional defocus value
  Real phaseShift(Real thetaX, Real thetaY, Real defocus = 0) const;

  /// Returns the transfer function @f$t=a * \exp(-\text{i}¸,k\,\chi)@f$ for the
  /// given scattering angles. In contrast to @ref pupilAngle() this function
  /// respects the aperture.
  Complex value(Real thetaX, Real thetaY, Real defocus = 0) const;

  /// Return the pupil function @f$p=\exp(-\text{i}¸,k\,\chi)@f$ for the given
  /// scattering angles. This function does not respect the aperture.
  Complex pupilAngle(Real thetaX, Real thetaY, Real defocus = 0) const
  {
    return fromArg(-phaseShift(thetaX, thetaY, defocus));
  }

  /// Returns the objective aperture in rad.
  Real aperture() const
  {
    return std::sqrt(mApertureSqr);
  }

  /// Returns the n-th astigmatism in Angstrom. The n-th astigmatism has an
  /// (n+1)-fold symmetry.
  Complex a(size_t order) const
  {
    return mA[order] * (order + R(1.0));
  }

  /// Returns the n-th order axial coma (B2, B4, B6) in Angstrom.
  Complex b(size_t order) const
  {
    assert(2 == order || 4 == order || 6 == order);
    return mB[(order - 2) / 2];
  }

  /// Returns the spherical aberrations (C1 = defocus, C3 = CS, C5, C7) in
  /// Angstrom.
  Real c(size_t order) const
  {
    assert(1 == order || 3 == order || 5 == order || 7 == order);
    switch (order)
    {
    case 1: return mC1;
    case 3: return mC3;
    case 5: return mSixthC5 * 6;
    case 7: return mEighthC7 * 8;
    default: return 0;
    }
  }

  /// Returns the n-th order three-lobe aberration (D4, D6) in Angstrom.
  Complex d(size_t order) const
  {
    assert(4 == order || 6 == order);
    if (4 == order)
      return mD4;
    else
      return mD6;
  }

  /// Returns the sixth-order pentacle aberration @f$F_6@f$ in Angstrom.
  Complex f(size_t order) const
  {
    (void)order;
    assert(6 == order);
    return mF6;
  }

  /// Returns the seventh-order chaplet aberration @f$G_7@f$ in Angstrom.
  Complex g(size_t order) const
  {
    (void)order;
    assert(7 == order);
    return mG7;
  }

  /// Returns the n-th order rosette aberration (R5, R7) in Angstrom.
  Complex r(size_t order) const
  {
    assert(5 == order || 7 == order);
    if (5 == order)
      return mR5;
    else
      return mR7;
  }

  /// Returns the n-th order star aberration (S3, S5, S7) in Angstrom.
  Complex s(size_t order) const
  {
    assert(3 == order || 5 == order || 7 == order);
    return mS[(order - 3) / 2];
  }

  Real scherzerDefocus(size_t n = 1) const;
  // todo: lichteDefocus
  // todo: lentzenDefocus

  /// Sets the objective aperture in rad.
  void setAperture(Real aperture)
  {
    mApertureSqr = aperture * aperture;
  }

  /// Sets the n-th astigmatism @f$A_n@f$ in Angstrom. The n-th astigmatism has
  /// a (n+1)-fold symmetry.
  void setA(size_t order, const Complex& a)
  {
    assert(order <= 7 && "Invalid order.");
    mA[order] = a / (order + R(1.0));
  }

  /// Sets the n-th order axial coma (B2, B4, B6) in Angstrom.
  void setB(size_t order, const Complex& value)
  {
    assert((2 == order || 4 == order || 6 == order) && "Invalid order.");
    mB[(order - 2) / 2] = value;
  }

  /// sets the spherical aberrations (C1 = defocus, C3 = CS, C5, C7) in
  /// Angstrom.
  void setC(size_t order, Real value)
  {
    assert(1 == order || 3 == order || 5 == order || 7 == order);
    switch (order)
    {
    case 1:
      mC1 = value;
      mHalfC1 = value / 2;
      break;
    case 3:
      mC3 = value;
      mFourthC3 = value / 4;
      break;
    case 5:
      mSixthC5 = value / 6;
      break;
    case 7:
      mEighthC7 = value / 8;
      break;
    }
  }

  /// Sets the n-th order three-lobe aberrations (D4, D6) in Angstrom.
  void setD(size_t order, const Complex& value)
  {
    assert(4 == order || 6 == order);
    if (4 == order)
      mD4 = value;
    else
      mD6 = value;
  }

  /// Sets the sixth-order pentacle aberration @f$F_6@f$ in Angstrom.
  void setF(size_t order, const Complex& f6)
  {
    (void)order;
    assert(6 == order);
    mF6 = f6;
  }

  /// Sets the seventh-order chaplet aberration @f$G_7@f$ in Angstrom.
  void setG(size_t order, const Complex& g7)
  {
    (void)order;
    assert(7 == order);
    mG7 = g7;
  }

  /// Sets the n-th order rosette aberration (R5, R7) in Angstrom.
  void setR(size_t order, const Complex& value)
  {
    assert(5 == order || 7 == order);
    if (5 == order)
      mR5 = value;
    else
      mR7 = value;
  }

  /// Sets the n-th order star aberration (S3, S5, S7) in Angstrom.
  void setS(size_t order, const Complex& value)
  {
    assert(3 == order || 5 == order || 7 == order);
    mS[(order - 3) / 2] = value;
  }

  /// This function returns the coherent CTF, i.e. the CTF without the
  /// damping envelopes.
  /// @param info
  ///  describes the properties of the target buffer.
  /// @param factor
  ///  Scaling factor for the whole CTF. This can be used to correct the FFT
  ///  scaling.
  /// @param defocus
  ///  Additional defocus value
  CBuffer2D buffer(const Buffer2DInfoBase& info, Real factor = 1,
                   Real defocus = 0) const;

  static std::shared_ptr<CtfCoherent> fromCommandLine(Real energy);
  void toCommandLine();

private:
  // the objective aperture in rad
  Real mApertureSqr = 1;

  // n-th astigmatism in Angstrom. These values are premultipled by 1 / (n + 1)
  // n-th astigmatism has an (n+1)-fold symmetry
  std::array<Complex, 8> mA = {{0, 0, 0, 0, 0, 0, 0, 0}};

  // axial coma (B2, B4, B6) in Angstrom
  std::array<Complex, 3> mB = {{0, 0, 0}};

  // half of the defocus in Angstrom
  Real mHalfC1 = 0;
  // defocus in Angstrom
  Real mC1 = 0;
  // one fourth of the third-order spherical aberration in Angstrom
  Real mFourthC3 = 0;
  // third-order spherical aberration in Angstrom
  Real mC3 = 0;
  // one sixth of the fifth-order spherical aberration in Angstrom
  Real mSixthC5 = 0;
  // one eighth of the seventh-order spherical aberration in Angstrom
  Real mEighthC7 = 0;

  // fourth-order three-lobe aberration in Angstrom
  Complex mD4 = 0;
  // sixth-order three-lobe aberration in Angstrom
  Complex mD6 = 0;

  // sixth-order pentacle aberration in Angstrom
  Complex mF6 = 0;

  // seventh-order chaplet aberration in Angstrom
  Complex mG7 = 0;

  // fifth-order rosette aberration in Angstrom
  Complex mR5 = 0;
  // seventh-order rosette aberration in Angstrom
  Complex mR7 = 0;

  // star aberration (S3, S5, S7) in Angstrom
  std::array<Complex, 3> mS = {{0, 0, 0}};
};

} // namespace Aurora

#endif
