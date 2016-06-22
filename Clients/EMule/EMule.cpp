//--- Aurora/Clients/Emule/EMule.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Bandlimit.hpp"
#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Ctf.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/PhaseShift.hpp"
#include "AuroraLib/Potential.hpp"
#include "AuroraLib/Propagator.hpp"
#include "AuroraLib/MultisliceParameters.hpp"
#include "AuroraLib/MultisliceOptions.hpp"
#include "AuroraLib/Progress.hpp"
#include "AuroraLib/Quadrature.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/StemDetector.hpp"
#include "AuroraLib/StemSimulation.hpp"
#include "AuroraLib/TemSimulation.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <iomanip>
#include <iostream>
#include <random>

namespace Aurora
{

enum class Mode
{
  tem, eftem, stem
};

// General
CL::Option<Mode> mode("mode", "", Mode::tem,
  enumVal(Mode::tem), enumVal(Mode::eftem), enumVal(Mode::stem));
CL::Option<std::string> outName("o", "Name of the output file", "NOTSET");
CL::Option<std::string> propagatorName("propagator", "", "Classical-FT");

// STEM
CL::Option<Real> alpha("alpha", "Inner ADF detector angle", R(0.05));
CL::Option<Real> beta("beta", "Outer ADF detector angle", R(0.15));

// Energy-filtered TEM (EFTEM)
CL::Option<Real> cc("CC", "C_C in Angstrom", R(1.5e7));
CL::Option<size_t> hermiteOrder("hermiteOrder", "", 10);
CL::Option<Real> energySpread("energySpread", "Energy spread in eV", R(0.5));
CL::Option<Real> highTensionRipple("highTensionRipple", "", R(1e-6));

CL::Option<bool> showProjPotential("showProjPotential", "", true);

// CBED
CL::Option<Real> convergenceAngle("convergenceAngle", "CBED convergence angle",
                                  0);

class EMule
{
public:
  EMule(int argc, char** argv);

private:
  void saveWave(const CBuffer2D& wave, const MultisliceParameters& params);
  void saveWave(const CBuffer2D& wave, const MultisliceParameters& params,
                size_t waveNum);

  void tem(const MultisliceParameters& params, const SamplePtr& sample);
  void eftem(MultisliceParameters params, const SamplePtr& sample);
  void stem(const MultisliceParameters& params, const SamplePtr& sample);

#if 0
  void onWave(int64_t index, const CBuffer2D& wave,
              const MultisliceParameters& params)
  {
    const int64_t reducedIndex = index % 2001;
    const int64_t iteration = index / 2001;

    if (1998 == reducedIndex)
      a = wave.clone();

    if (1999 == reducedIndex)
      b = wave.clone();

    if (2000 == reducedIndex)
    {
      c = wave.clone();

      CBuffer2D firstDerivative(wave.cols(), wave.rows());
      firstDerivative.assign(R(2.0) * params.wavenumber() / params.deltaZ() * (c - b));
      CBuffer2D secondDerivative(wave.cols(), wave.rows());
      secondDerivative.assign(R(1.0) / powerOf<2>(params.deltaZ()) * (c - R(2.0) * b + a));

      firstDerivative.save(params.debugDir() + "/1st-" + std::to_string(iteration) + ".tiff", params, "");
      secondDerivative.save(params.debugDir() + "/2nd-" + std::to_string(iteration) + ".tiff", params, "");
    }
  }
#endif

  CBuffer2D a;
  CBuffer2D b;
  CBuffer2D c;
};

EMule::EMule(int argc, char** argv)
{
  std::cerr << "EMule Client\n";

  if (verbosity >= 1)
    std::cout << buildString() << '\n';

  // Add the form factor names to the description text
  formFactorName.setDescription(FormFactorParametrization::namesString());
  phaseShiftName.setDescription(PhaseShift::namesString());
  potentialName.setDescription(Potential::namesString());
  propagatorName.setDescription(Propagator::namesString());

  // parse the command line
  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
  {
    std::cout << "Simulation started\n";

    // output the command line options
    std::cout << CL::configString() << '\n';
  }

  MultisliceParameters params = multisliceParamsFromCommandLine();

  // load the sample description
  auto sample = std::make_shared<Sample>(sampleName);

  for (auto& atom : sample->atoms())
    atom.setBFactor(bFactor);

  if (Mode::tem == mode)
    tem(params, sample);
  else if (Mode::eftem == mode)
    eftem(params, sample);
  else if (Mode::stem == mode)
    stem(params, sample);
}

void EMule::saveWave(const CBuffer2D& wave, const MultisliceParameters& params)
{
  std::ostringstream description;
  description << "# EMule, Version: " << version << "\n"
                 "# " << Time::now() << '\n'
              << CL::configString() << '\n';

  wave.save(outName, params, description.str());
}

void EMule::saveWave(const CBuffer2D& wave, const MultisliceParameters& params,
                     size_t waveNum)
{
  std::ostringstream description;
  description << "# EMule, Version: " << version << "\n"
                 "# " << Time::now() << '\n'
              << CL::configString() << '\n';

  wave.save(outName + std::to_string(waveNum), params, description.str());
}

void EMule::tem(const MultisliceParameters& params, const SamplePtr& sample)
{
  if (verbosity >= 1)
    std::cout << "TEM mode\n";

  TemSimulation tem(params, numRepetitions, sample, propagatorName);

  using namespace std::placeholders;
  if (isTerminal(std::cout))
    tem.onWave = std::bind(&progress, _1, _2);

  // activate the next line, if you want to output the first and second
  // derivatives of the wave
#if 0
  tem.onWave = std::bind(&EMule::onWave, this, _1, _3, params);
#endif

  if (0 != convergenceAngle)
  {
    std::cout << "CBED\n";
    tem.setCbed(convergenceAngle);
  }

  // execute the simulation
  tem.run();

  saveWave(tem.wave(), params);
}

void EMule::eftem(MultisliceParameters params, const SamplePtr& sample)
{
  if (verbosity >= 1)
    std::cout << "EFTEM mode\n";

  HermitePolynomial polynomial(hermiteOrder);
  Real focusSpread = cc * std::sqrt(powerOf<2>(energySpread / energy) +
                                    powerOf<2>(highTensionRipple));
  std::cout << "Focus spread = " << focusSpread * R(0.1) << " nm";

  for (size_t i = 0; hermiteOrder != i; ++i)
  {
    Real focusShift = polynomial.root(i) * Math::sqrt2 * focusSpread;
    Real energyShift = focusShift / cc * energy;

    params.setEnergy(energy + energyShift);

    TemSimulation tem(params, numRepetitions, sample, propagatorName);

    // execute the simulation
    tem.run();

    saveWave(tem.wave(), params, i);
  }
}

void EMule::stem(const MultisliceParameters& params, const SamplePtr& sample)
{
  if (verbosity >= 1)
    std::cout << "STEM mode\n";

  auto ctf = CtfCoherent::fromCommandLine(energy);
  auto detector = std::make_shared<AdfDetector>(params, alpha, beta);

  StemSimulation stem(params, numRepetitions, sample, ctf, propagatorName,
                      detector);

  // execute the simulation
  stem.run();

  saveWave(stem.result(), params);
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::EMule emule(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
