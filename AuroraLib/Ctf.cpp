//--- Aurora/AuroraLib/Ctf.cpp -------------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Ctf.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Physics.hpp"

namespace Aurora
{

// coherent CTF parameters
namespace CtfParameters
{
//ctf.setAperture(options.as<Real>("aperture") * R(1e3));
CL::Option<Real> aperture("aperture", "Objective aperture in rad", 1.0);

CL::Option<Real> a0Abs("A0Abs", "Modulus of the beam/image shift in Angstrom");
CL::Option<Real> a0Ang("A0Ang", "Angle of the beam/image shift in degree");
CL::Option<Real> a1Abs("A1Abs", "Modulus of the twofold astigmatism in Angstrom");
CL::Option<Real> a1Ang("A1Ang", "Angle of the twofold astigmatism in degree");
CL::Option<Real> a2Abs("A2Abs", "Modulus of the threefold astigmatism in Angstrom");
CL::Option<Real> a2Ang("A2Ang", "Angle of the threefold astigmatism in degree");
CL::Option<Real> a3Abs("A3Abs", "Modulus of the fourfold astigmatism in Angstrom");
CL::Option<Real> a3Ang("A3Ang", "Angle of the fourfold astigmatism in degree");
CL::Option<Real> a4Abs("A4Abs", "Modulus of the fivefold astigmatism in Angstrom");
CL::Option<Real> a4Ang("A4Ang", "Angle of the fivefold astigmatism in degree");
CL::Option<Real> a5Abs("A5Abs", "Modulus of the sixfold astigmatism in Angstrom");
CL::Option<Real> a5Ang("A5Ang", "Angle of the sixfold astigmatism in degree");
CL::Option<Real> a6Abs("A6Abs", "Modulus of the sevenfold astigmatism in Angstrom");
CL::Option<Real> a6Ang("A6Ang", "Angle of the sevenfold astigmatism in degree");
CL::Option<Real> a7Abs("A7Abs", "Modulus of the eightfold astigmatism in Angstrom");
CL::Option<Real> a7Ang("A7Ang", "Angle of the eightfold astigmatism in degree");

CL::Option<Real> b2Abs("B2Abs", "Modulus of the second-order axial coma in Angstrom");
CL::Option<Real> b2Ang("B2Ang", "Angle of the second-order axial coma in degree");
CL::Option<Real> b4Abs("B4Abs", "Modulus of the fourth-order axial coma in Angstrom");
CL::Option<Real> b4Ang("B4Ang", "Angle of the fourth-order axial coma in degree");
CL::Option<Real> b6Abs("B6Abs", "Modulus of the sixth-order axial coma in Angstrom");
CL::Option<Real> b6Ang("B6Ang", "Angle of the sixth-order axial coma in degree");

CL::Option<Real> c1("C1", "Defocus in Angstrom");
CL::Option<Real> c3("C3", "Third-order spherical aberration in Angstrom");
CL::Option<Real> c5("C5", "Fifth-order spherical aberration in Angstrom");
CL::Option<Real> c7("C7", "Seventh-order spherical aberration in Angstrom");

CL::Option<Real> d4Abs("D4Abs", "Modulus of the fourth-order three-lobe aberration in Angstrom");
CL::Option<Real> d4Ang("D4Ang", "Angle of the fourth-order three-lobe aberration in degree");
CL::Option<Real> d6Abs("D6Abs", "Modulus of the sixth-order three-lobe aberration in Angstrom");
CL::Option<Real> d6Ang("D6Ang", "Angle of the sixth-order three-lobe aberration in degree");

CL::Option<Real> f6Abs("F6Abs", "Modulus of the sixth-order pentacle aberration aberration in Angstrom");
CL::Option<Real> f6Ang("F6Ang", "Angle of the sixth-order pentacle aberration aberration in degree");

CL::Option<Real> g7Abs("G7Abs", "Modulus of the seventh-order chaplet aberration aberration in Angstrom");
CL::Option<Real> g7Ang("G7Ang", "Angle of the seventh-order chaplet aberration aberration in degree");

CL::Option<Real> r5Abs("R5Abs", "Modulus of the fifth-order rosette aberration aberration in Angstrom");
CL::Option<Real> r5Ang("R5Ang", "Angle of the fifth-order rosette aberration aberration in degree");
CL::Option<Real> r7Abs("R7Abs", "Modulus of the seventh-order rosette aberration aberration in Angstrom");
CL::Option<Real> r7Ang("R7Ang", "Angle of the seventh-order rosette aberration aberration in degree");

CL::Option<Real> s3Abs("S3Abs", "Modulus of the third-order star aberration aberration in Angstrom");
CL::Option<Real> s3Ang("S3Ang", "Angle of the third-order star aberration aberration in degree");
CL::Option<Real> s5Abs("S5Abs", "Modulus of the fifth-order star aberration aberration in Angstrom");
CL::Option<Real> s5Ang("S5Ang", "Angle of the fifth-order star aberration aberration in degree");
CL::Option<Real> s7Abs("S7Abs", "Modulus of the seventh-order star aberration aberration in Angstrom");
CL::Option<Real> s7Ang("S7Ang", "Angle of the seventh-order star aberration aberration in degree");
//  ctf.setEnergy(options.as<Real>("kineticEnergy") * R(1e3));

} // namespace CtfParameters

Real CtfCoherent::scherzerDefocus(size_t n) const
{
  return -std::copysign(std::sqrt((2 * n - R(0.5)) * std::abs(mC3) *
                                  wavelength()), mC3);
}

Complex CtfCoherent::value(Real thetaX, Real thetaY, Real defocus) const
{
  // the squared length of the theta-vector
  const Real thetaSqr = thetaX * thetaX + thetaY * thetaY;

  if (thetaSqr > mApertureSqr)
    return Complex(0.0);

  return fromArg(-phaseShift(thetaX, thetaY, defocus));
}

Real CtfCoherent::phaseShift(Real thetaX, Real thetaY, Real defocus) const
{
  // the complex conjugated angle omega
  const Real omegaNorm = thetaX * thetaX + thetaY * thetaY;
  const Complex omegaConj(thetaX, -thetaY);
  const Complex omegaConj2 = omegaConj  * omegaConj;
  const Complex omegaConj3 = omegaConj2 * omegaConj;
  const Complex omegaConj4 = omegaConj2 * omegaConj2;
  const Complex omegaConj5 = omegaConj3 * omegaConj2;
  const Complex omegaConj6 = omegaConj3 * omegaConj3;
  const Complex omegaConj7 = omegaConj4 * omegaConj3;

  Real chi = (omegaConj * (mA[0] + mA[1] * omegaConj + mB[0] * conj(omegaConj2)
    + mA[2] * omegaConj2 + mS[0] * conj(omegaConj3) + mA[3] * omegaConj3
    + mD4 * conj(omegaConj4) + mB[1] * conj(omegaConj3) * omegaConj
    + mA[4] * omegaConj4 + mR5 * conj(omegaConj5)
    + mS[1] * conj(omegaConj4) * omegaConj + mA[5] * omegaConj5
    + mF6 * conj(omegaConj6) + mD6 * conj(omegaConj5) * omegaConj
    + mB[2] * conj(omegaConj4) * omegaConj2 + mA[6] * omegaConj6
    + mG7 * conj(omegaConj7) + mR7 * conj(omegaConj6) * omegaConj
    + mS[2] * conj(omegaConj5) * omegaConj2 + mA[7] * omegaConj7)).real();

  chi += omegaNorm * (mHalfC1 + R(0.5) * defocus +
         omegaNorm * (mFourthC3 +
         omegaNorm * (mSixthC5 +
         omegaNorm * mEighthC7)));

  return chi * wavenumber();
}

CBuffer2D CtfCoherent::buffer(const Buffer2DInfoBase& info, Real factor,
                      Real defocus) const
{
  const int64_t cols = info.cols();
  const int64_t rows = info.rows();

  CBuffer2D result(cols, rows);

  const Real deltaThetaX = info.deltaKX() / wavenumber();
  const Real deltaThetaY = info.deltaKY() / wavenumber();

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
  {
    // y coordinate in Fourier space
    const Real thetaY = deltaThetaY * (y > rows / 2 ? y - rows : y);

    for (int64_t x = 0; x < cols; ++x)
    {
      // x coordinate in Fourier space
      const Real thetaX = deltaThetaX * (x > cols / 2 ? x - cols : x);
      result.pixel(x, y) = factor * value(thetaX, thetaY, defocus);
    }
  }

  return result;
}

/* static */ std::shared_ptr<CtfCoherent> CtfCoherent::fromCommandLine(Real energy)
{
  auto result = std::make_shared<CtfCoherent>(energy);

  using namespace CtfParameters;
  result->setAperture(CtfParameters::aperture);

  result->setA(0, std::polar<Real>(a0Abs, a0Ang * Math::degToRad));
  result->setA(1, std::polar<Real>(a1Abs, a1Ang * Math::degToRad));
  result->setA(2, std::polar<Real>(a2Abs, a2Ang * Math::degToRad));
  result->setA(3, std::polar<Real>(a3Abs, a3Ang * Math::degToRad));
  result->setA(4, std::polar<Real>(a4Abs, a4Ang * Math::degToRad));
  result->setA(5, std::polar<Real>(a5Abs, a5Ang * Math::degToRad));
  result->setA(6, std::polar<Real>(a6Abs, a6Ang * Math::degToRad));
  result->setA(7, std::polar<Real>(a7Abs, a7Ang * Math::degToRad));

  result->setB(2, std::polar<Real>(b2Abs, b2Ang * Math::degToRad));
  result->setB(4, std::polar<Real>(b4Abs, b4Ang * Math::degToRad));
  result->setB(6, std::polar<Real>(b6Abs, b6Ang * Math::degToRad));

  result->setC(1, CtfParameters::c1);
  result->setC(3, CtfParameters::c3);
  result->setC(5, CtfParameters::c5);
  result->setC(7, CtfParameters::c7);

  result->setD(4, std::polar<Real>(d4Abs, d4Ang * Math::degToRad));
  result->setD(6, std::polar<Real>(d6Abs, d6Ang * Math::degToRad));

  result->setF(6, std::polar<Real>(f6Abs, f6Ang * Math::degToRad));

  result->setG(7, std::polar<Real>(g7Abs, g7Ang * Math::degToRad));

  result->setR(5, std::polar<Real>(r5Abs, r5Ang * Math::degToRad));
  result->setR(7, std::polar<Real>(r7Abs, r7Ang * Math::degToRad));

  result->setS(3, std::polar<Real>(s3Abs, s3Ang * Math::degToRad));
  result->setS(5, std::polar<Real>(s5Abs, s5Ang * Math::degToRad));
  result->setS(7, std::polar<Real>(s7Abs, s7Ang * Math::degToRad));

  return result;
}

void CtfCoherent::toCommandLine()
{
  using namespace CtfParameters;

  CtfParameters::aperture.setValue(aperture());

  a0Abs.setValue(std::abs(a(0)));
  a0Ang.setValue(std::arg(a(0)));
  a1Abs.setValue(std::abs(a(1)));
  a1Ang.setValue(std::arg(a(1)));
  a2Abs.setValue(std::abs(a(2)));
  a2Ang.setValue(std::arg(a(2)));
  a3Abs.setValue(std::abs(a(3)));
  a3Ang.setValue(std::arg(a(3)));
  a4Abs.setValue(std::abs(a(4)));
  a4Ang.setValue(std::arg(a(4)));
  a5Abs.setValue(std::abs(a(5)));
  a5Ang.setValue(std::arg(a(5)));
  a6Abs.setValue(std::abs(a(6)));
  a6Ang.setValue(std::arg(a(6)));
  a7Abs.setValue(std::abs(a(7)));
  a7Ang.setValue(std::arg(a(7)));

  b2Abs.setValue(std::abs(b(2)));
  b2Ang.setValue(std::arg(b(2)));
  b4Abs.setValue(std::abs(b(4)));
  b4Ang.setValue(std::arg(b(4)));
  b6Abs.setValue(std::abs(b(6)));
  b6Ang.setValue(std::arg(b(6)));

  CtfParameters::c1.setValue(c(1));
  CtfParameters::c3.setValue(c(3));
  CtfParameters::c5.setValue(c(5));
  CtfParameters::c7.setValue(c(7));

  d4Abs.setValue(std::abs(d(4)));
  d4Ang.setValue(std::arg(d(4)));
  d6Abs.setValue(std::abs(d(6)));
  d6Ang.setValue(std::arg(d(6)));

  f6Abs.setValue(std::abs(f(6)));
  f6Ang.setValue(std::arg(f(6)));

  g7Abs.setValue(std::abs(g(7)));
  g7Ang.setValue(std::arg(g(7)));

  r5Abs.setValue(std::abs(r(5)));
  r5Ang.setValue(std::arg(r(5)));
  r7Abs.setValue(std::abs(r(7)));
  r7Ang.setValue(std::arg(r(7)));

  s3Abs.setValue(std::abs(s(3)));
  s3Ang.setValue(std::arg(s(3)));
  s5Abs.setValue(std::abs(s(5)));
  s5Ang.setValue(std::arg(s(5)));
  s7Abs.setValue(std::abs(s(7)));
  s7Ang.setValue(std::arg(s(7)));
}

} // namespace Aurora
