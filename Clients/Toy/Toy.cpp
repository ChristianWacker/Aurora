//--- Aurora/Clients/Toy/Toy.cpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Physics.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <Eigen/Core>
#include <fstream>
#include <iostream>

namespace Aurora
{

typedef Eigen::Matrix<Complex, 2, 2> Matrix2;

std::ofstream phaseShiftFile;
std::ofstream b21File;

enum class Mode
{
  rect, triangle, gaussian, peng, cosh
};

CL::Option<Mode> mode("mode", "", Mode::gaussian,
  enumVal(Mode::rect), enumVal(Mode::triangle), enumVal(Mode::gaussian),
  enumVal(Mode::peng), enumVal(Mode::cosh));
CL::Option<Real> energy("energy", "Energy of the electron in eV", 2e4);
CL::Option<Real> v0("v0", "", 1.0);
CL::Option<Real> a("a", "Width of the potential in Angstrom", 1.0);
CL::Option<int> maxNumSlices("maxNumSlices", "", 1500);

/// Class to calculate the average value of a potential on a certain interval
/// [x1, x2]
class AveragedPotential
{
public:
  virtual ~AveragedPotential() {}
  /// The antiderivative of the potential
  virtual Real antiderivative(Real) const = 0;

  Real eval(Real x1, Real x2) const
  {
    return v0 * (antiderivative(x2) - antiderivative(x1)) / (x2 - x1);
  }
};

/// Square with width 2 a and height 1
class PotentialRect : public AveragedPotential
{
public:
  Real antiderivative(Real x) const override
  {
    if (x < -1 / a)
      return 0;
    else if (x < 1 / a)
      return -x - 1 / a;
    else
      return -2 / a;
  }
};

/// Triangle with width 1 and height 1
class PotentialTriangle : public AveragedPotential
{
public:
  Real antiderivative(Real x) const override
  {
    if (x < -R(0.5))
      return 0;
    else if (x < 0)
      return -x * x - x - R(0.25);
    else if (x < R(0.5))
      return x * x - x - R(0.25);
    else
      return -0.5;
  }
};

/// Gaussian function with standard deviation a
class PotentialGaussian : public AveragedPotential
{
public:
  Real antiderivative(Real x) const override
  {
    return -(1 + std::erf(x / a / Math::sqrt2));
  }
};

/// Potential: -1 / Cosh(a * x)^2
class PotentialCosh : public AveragedPotential
{
public:
  Real antiderivative(Real x) const override
  {
    return -std::tanh(a * x) / a;
  }
};

/// Potential based on the form factor parametrization of peng. As the form
/// factor parametrization is based on Guassian function it is safe to
/// evaluate the phase shift for very small radii.
class PotentialPeng : public AveragedPotential
{
public:
  PotentialPeng() :
    mTable(FormFactorParametrization::create("Peng")),
    mElement(mTable->element(14, 0))
  {}

  Real antiderivative(Real x) const override
  {
    using namespace Phys;
    return -hBar * hBar / mE * mElement.phaseShift(0.0, x, 0.0);
  }

private:
  std::unique_ptr<FormFactorParametrization> mTable;
  const Element& mElement;
};

/// Application class
class Toy
{
public:
  Toy(int argc, char** argv);

private:
  /// Return a transfer matrix transforming the coefficients for the wave
  /// function on the interval [x1, x2] to the coefficients for the interval
  /// [x2, x3]
  Matrix2 matrixForInterval(Real x1, Real x2, Real x3);

  /// Return a transfer matrix for the whole interval based on numSlices slices.
  Matrix2 transferMatrix(int numSlices);

  std::unique_ptr<AveragedPotential> mPotential;
};

Toy::Toy(int argc, char** argv)
{
  std::cerr << "Toy Client\n";

  // parse the command line
  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
  {
    std::cout << buildString() << '\n';

    // output the command line options
    std::cout << CL::configString() << '\n';
  }

  if (Mode::rect == mode)
    mPotential = std::make_unique<PotentialRect>();
  else if (Mode::triangle == mode)
    mPotential = std::make_unique<PotentialTriangle>();
  else if (Mode::gaussian == mode)
    mPotential = std::make_unique<PotentialGaussian>();
  else if (Mode::peng == mode)
    mPotential = std::make_unique<PotentialPeng>();
  else if (Mode::cosh == mode)
    mPotential = std::make_unique<PotentialCosh>();

  for (int numSlices = 2; numSlices < maxNumSlices; ++numSlices)
  {
    std::string postfix = std::to_string(numSlices) + "slices.txt";
    phaseShiftFile.open("phaseShift-" + postfix);
    b21File.open("b21-" + postfix);
    auto mat = transferMatrix(numSlices);
    std::cout << numSlices << ' ' << absSqr(mat(1, 0) / mat(1, 1)) << '\n';
 // std::cout << numSlices << ' ' << abs(transferMatrix(numSlices)(1, 0)) << '\n';
    b21File.close();
    phaseShiftFile.close();
  }
}

Matrix2 Toy::matrixForInterval(Real x1, Real x2, Real x3)
{
  using namespace Phys;

  // slice thickness
  const Real deltaZ = x2 - x1;

  phaseShiftFile << mPotential->eval(x1, x2) << '\n';

  constexpr Real factor = 2 * mE / powerOf<2>(hBar);
  // Wavenumbers on the two intervals
  const Real k1 = std::sqrt(factor * (energy - mPotential->eval(x1, x2)));
  const Real k2 = std::sqrt(factor * (energy - mPotential->eval(x2, x3)));

  Matrix2 mat;
  mat(0, 0) = R(0.5) * (1 + k2 / k1) * exp(Complex(0,  k1 * deltaZ));
  mat(0, 1) = R(0.5) * (1 - k2 / k1) * exp(Complex(0, -k1 * deltaZ));
  mat(1, 0) = R(0.5) * (1 - k2 / k1) * exp(Complex(0,  k1 * deltaZ));
  mat(1, 1) = R(0.5) * (1 + k2 / k1) * exp(Complex(0, -k1 * deltaZ));

  return mat;
}

Matrix2 Toy::transferMatrix(int numSlices)
{
  assert(numSlices >= 2);

  const Real start = -1.5;
  const Real end = 1.5;
  const Real deltaZ = (end - start) / numSlices;

  // start with the identity matrix
  Matrix2 result = Matrix2::Identity();

  // for every slices
  for (int i = -1; i < numSlices; ++i)
  {
    const Real x1 = start + i * deltaZ;
    const Real x2 = start + (i + 1) * deltaZ;
    const Real x3 = start + (i + 2) * deltaZ;
    Matrix2 mat = matrixForInterval(x1, x2, x3);
    // combine the previous transfer matrix with the new one
    result = mat * result;

    b21File << std::abs(result(1, 0)) << '\n';
  }

  return result;
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Toy toy(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
