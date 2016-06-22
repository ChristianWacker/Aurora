//--- Aurora/AuroraLib/MicroscopeSimulation.cpp --------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/MicroscopeSimulation.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/Version.hpp"

#include <iostream>
#include <random>

namespace Aurora
{

CL::Option<int> seed("seed", "seed for the random number generator", -1);

ConstSamplePtr MicroscopeSimulation::displacedSample() const
{
  assert(mParams.thermicScattering() == ThermicScattering::phonon);

  // initialize the random number generator
  std::random_device randomDevice;

  std::mt19937 generator;
  if (-1 == seed)
    generator.seed(randomDevice());
  else
    generator.seed(seed.value());

  SamplePtr result = std::make_shared<Sample>();

  // move the atoms around
  for (const auto& atom : mSample->atoms())
  {
    const Real meanDisplacement = std::sqrt(2 * atom.bFactor());
    std::normal_distribution<Real> distribution(0, meanDisplacement);

    const Vector3R displacementVector(distribution(generator),
                                      distribution(generator),
                                      distribution(generator));

    result->atoms().push_back(Atom(atom.position() + displacementVector,
                                   atom.atomicNumber(), atom.charge()));
  }

  result->preprocess();
  return result;
}

void MicroscopeSimulation::printStartupMessages() const
{
  if (verbosity >= 2)
  {
    // output the dimensions of simulated volume
    std::cout << "Simulation Volume:\n"
                 "minX = " << mParams.minBox().x() << "\n"
                 "maxX = " << mParams.maxBox().x() << "\n"
                 "minY = " << mParams.minBox().y() << "\n"
                 "maxY = " << mParams.maxBox().y() << "\n"
                 "minZ = " << mParams.minBox().z() << "\n"
                 "maxZ = " << mParams.maxBox().z() << '\n';
  }

  if (verbosity >= 1)
  {
    std::cout << "Wavelength (lambda) = "
              << mParams.wavelength() * R(100.0) << " pm\n"
                 "Relativistic factor (gamma) = "
              << mParams.lorentzFactor() << "\n"
                 "Speed = "
              << mParams.speed() << " c\n"
                 "Nyquist frequency in x-Direction = "
              << mParams.kXNyquist() / Math::twoPi << "1 / Angstrom\n"
                 "Nyquist frequency in y-Direction = "
              << mParams.kYNyquist() / Math::twoPi << "1 / Angstrom"
              << std::endl;
  }
}

} // namespace Aurora
