//--- Aurora/Clients/Matrix/Matrix.cpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "Matrix.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/MultisliceOptions.hpp"
#include "AuroraLib/PhaseShift.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/StemDetector.hpp"
#include "AuroraLib/StemSimulation.hpp"
#include "AuroraLib/Version.hpp"

#include "Pass1.hpp"
#include "Pass2.hpp"
#include "Pass3.hpp"
#include "PassBandlimit.hpp"
#include "PassPower.hpp"
#include "PassScanning.hpp"
#include "PassTem.hpp"

#include <iomanip>
#include <iostream>

namespace Aurora
{

enum class Mode
{
  laplace, matrixA2, pass1, pass2F, pass2FLarge, pass2T, pass2TLarge,
  pass2TSLarge, pass3F, pass3T, pass3TLarge, pass3SLarge, powerF, powerT,
  powerTLarge, powerSLarge, iterativeTem, tem, temLarge, scanning,
  iterativeScanning, bandlimit, simulateScanning, inspect
};

CL::Option<std::string> inTemplate("i", "Name template for the input files", "NOTSET");
CL::Option<std::string> outName("o", "Name of the output file", "NOTSET");
CL::Option<size_t> part("part", "", 1);
CL::Option<size_t> numParts("numParts", "", 1);
CL::Option<Mode> mode("mode", "", Mode::pass1,
  enumVal(Mode::laplace), enumVal(Mode::matrixA2), enumVal(Mode::pass1),
  enumVal(Mode::pass2F), enumVal(Mode::pass2FLarge), enumVal(Mode::pass2T),
  enumVal(Mode::pass2TLarge), enumVal(Mode::pass2TSLarge),
  enumVal(Mode::pass3F), enumVal(Mode::pass3T), enumVal(Mode::pass3TLarge),
  enumVal(Mode::pass3SLarge), enumVal(Mode::powerF), enumVal(Mode::powerT),
  enumVal(Mode::powerTLarge), enumVal(Mode::powerSLarge), enumVal(Mode::tem),
  enumVal(Mode::temLarge), enumVal(Mode::iterativeTem), enumVal(Mode::scanning),
  enumVal(Mode::iterativeScanning), enumVal(Mode::bandlimit),
  enumVal(Mode::simulateScanning), enumVal(Mode::inspect));
CL::Option<std::string> matrixDir("matrixDir", "Directory for the matrices", "NOTSET");
CL::Option<bool> cacheA2("cacheA2", "", false);

// Power
CL::Option<size_t> increment("increment", "power mode", 1);
CL::Option<size_t> power("power", "power mode", 1);
CL::Option<size_t> startIndex("startIndex", "power mode", 1);

// Inspect
CL::Option<int> submatrix("submatrix", "", 0);

// SEM / STEM
CL::Option<Real> alpha("alpha", "Inner ADF detector angle", R(0.05));
CL::Option<Real> beta("beta", "Outer ADF detector angle", R(0.15));
CL::Option<bool> signals("outputSignals", "Write the detector signals in scanning mode into files", false);
CL::Option<int64_t> scanCols("scanCols", "", 128);
CL::Option<int64_t> scanRows("scanRows", "", 128);
CL::Option<int64_t> minIterations("minIterations", "", 9);
CL::Option<int64_t> maxIterations("maxIterations", "", 32);

Application::Application(int argc, char** argv)
{
  std::cerr << "Matrix Client\n";

  // parse the command line
  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
  {
    std::cout << buildString() << '\n';

    // output the command line options
    std::cout << CL::configString() << '\n' << std::endl;
  }

  // set the simulation parameters
  MultisliceParameters params = multisliceParamsFromCommandLine();

//  mParams.setDebug(Debug::pass2);

  if (verbosity >= 1)
  {
    std::cout << "Wavelength (lambda) = "
              << params.wavelength() * R(100.0) << " pm\n"
                 "Relativistic factor (gamma) = "
              << params.lorentzFactor() << "\n"
                 "Speed = "
              << params.speed() << " c\n"
                 "Nyquist frequency in x-Direction = "
              << params.kXNyquist() / Math::twoPi << "1 / Angstrom\n"
                 "Nyquist frequency in y-Direction = "
              << params.kYNyquist() / Math::twoPi << "1 / Angstrom"
              << std::endl;
  }

  if (numParts < 1)
    AURORA_THROW(ECommandLineError, "Number of parts must be one or greater");

  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  SamplePtr sample;
  // load the sample if necessary
  if ((Mode::matrixA2 == mode) || (Mode::pass1 == mode) ||
      (!cacheA2 && ((Mode::pass2F == mode) || (Mode::pass2FLarge == mode) ||
                    (Mode::pass2T == mode) || (Mode::pass2TLarge == mode) ||
                    (Mode::pass2TSLarge == mode))))
  {
    sample = std::make_shared<Sample>(sampleName);

    for (auto& atom : sample->atoms())
      atom.setBFactor(bFactor);
  }

  std::unique_ptr<Pass> pass;

  // numParts for pass2 and pass3 need to be the same. numParts for pass1 can
  // be different.
  if (Mode::laplace == mode)
  {
    pass = std::make_unique<PassLaplace>(params, matrixDir);
  }
  else if (Mode::matrixA2 == mode)
  {
    pass = std::make_unique<PassMatrixA2>(params, matrixDir, sample, part);
  }
  else if (Mode::pass1 == mode)
  {
    pass = std::make_unique<Pass1>(params, matrixDir, sample, numParts, part);
  }
  else if (Mode::pass2F == mode)
  {
    pass = std::make_unique<Pass2F>(params, matrixDir, sample, numParts, part);
  }
  else if (Mode::pass2FLarge == mode)
  {
    pass = std::make_unique<Pass2FLarge>(params, matrixDir, sample, numParts,
                                         part);
  }
  else if (Mode::pass2T == mode)
  {
    pass = std::make_unique<Pass2T>(params, matrixDir, sample, numParts, part);
  }
  else if (Mode::pass2TLarge == mode)
  {
    pass = std::make_unique<Pass2TLarge>(params, matrixDir, sample, numParts,
                                         part);
  }
  else if (Mode::pass2TSLarge == mode)
  {
    if (Bandlimit::Disabled == params.bandlimit())
    {
      pass = std::make_unique<Pass2TSLarge>(params, matrixDir, sample, numParts,
                                            part);
    }
    else
    {
      pass = std::make_unique<Pass2TSBandlimitLarge>(params, matrixDir, sample,
                                                     numParts, part);
    }
  }
  else if (Mode::pass3F == mode)
  {
    pass = std::make_unique<Pass3>(params, matrixDir, numParts, "F",
                                   Direction::forward);
  }
  else if (Mode::pass3T == mode)
  {
    pass = std::make_unique<Pass3>(params, matrixDir, numParts, "T",
                                   Direction::forward);
  }
  else if (Mode::pass3TLarge == mode)
  {
    pass = std::make_unique<Pass3Large>(params, matrixDir, numParts, "T",
                                        Direction::forward);
  }
  else if (Mode::pass3SLarge == mode)
  {
    pass = std::make_unique<Pass3Large>(params, matrixDir, numParts, "S",
                                        Direction::backward);
  }
  else if (Mode::powerF == mode)
  {
    pass = std::make_unique<PassPowerF>(params, matrixDir, power,
                                        outName);
  }
  else if (Mode::powerT == mode)
  {
    pass = std::make_unique<PassPowerT>(params, matrixDir, power);
  }
  else if (Mode::powerTLarge == mode)
  {
    pass = std::make_unique<PassPowerLarge>(params, matrixDir, power, "T",
                                            startIndex, increment);
  }
  else if (Mode::powerSLarge == mode)
  {
    pass = std::make_unique<PassPowerLarge>(params, matrixDir, power, "S",
                                            startIndex, increment);
  }
  else if (Mode::tem == mode)
  {
    pass = std::make_unique<PassTem>(params, matrixDir, power, outName);
  }
  else if (Mode::temLarge == mode)
  {
    pass = std::make_unique<PassTemLarge>(params, matrixDir, power, outName);
  }
  else if (Mode::iterativeTem == mode)
  {
    pass = std::make_unique<PassTemIterative>(params, matrixDir, power,
      outName, minIterations, maxIterations);
  }
  else if (Mode::scanning == mode)
  {
    auto stemDetector = std::make_shared<AdfDetector>(params, alpha, beta);
    auto passScanning = std::make_unique<PassScanning>(params, matrixDir,
      power, outName, stemDetector);

    if (signals)
    {
//      passScanning->setSignals(Signal::all);
      passScanning->setSignals(Signal::sem);
      passScanning->setSignals(Signal::stem);
      passScanning->setSignals(Signal::psiEx);
      passScanning->setSignals(Signal::psiBack);
    }
    else
    {
      passScanning->setSignals(Signal::sem);
      passScanning->setSignals(Signal::stem);
    }

    pass = std::move(passScanning);
  }
  else if (Mode::iterativeScanning == mode)
  {
    auto stemDetector = std::make_shared<AdfDetector>(params, alpha, beta);
    auto passScanning = std::make_unique<PassScanningIterative>(params,
      matrixDir, power, outName, stemDetector, scanCols, scanRows,
      minIterations, maxIterations);

    if (signals)
    {
  //    pass->setSignals(Signal::all);
      passScanning->setSignals(Signal::sem);
      passScanning->setSignals(Signal::stem);
      passScanning->setSignals(Signal::psiEx);
      passScanning->setSignals(Signal::psiBack);
    }
    else
    {
      passScanning->setSignals(Signal::sem);
      passScanning->setSignals(Signal::stem);
    }

    pass = std::move(passScanning);
  }
  else if (Mode::bandlimit == mode)
  {
    pass = std::make_unique<PassBandlimit>(params, matrixDir);
  }

  if (pass)
    pass->run();
  else if (Mode::simulateScanning == mode)
    simulateStem(inTemplate, outName); // todo: move this mode to another client
  else if (Mode::inspect == mode)
    inspect(inTemplate, submatrix);
  else
    AURORA_UNREACHABLE;

  if (verbosity >= 1)
  {
    std::cout  << "Total time needed: "
               << duration(high_resolution_clock::now() - startTime)
               << std::endl;
  }
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Application application(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
