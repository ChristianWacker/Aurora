//--- Aurora/Clients/Ice/Ice.cpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Bandlimit.hpp"
#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/MultisliceParameters.hpp"
#include "AuroraLib/MultisliceOptions.hpp"
#include "AuroraLib/Progress.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/TemSimulation.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include "InelasticSimulation.hpp"

#include <iomanip>
#include <iostream>
#include <random>

namespace Aurora
{

enum class Mode
{
  tem
};

// General
CL::Option<Mode> mode("mode", "", Mode::tem, enumVal(Mode::tem));
CL::Option<std::string> outName("o", "Name of the output file", "NOTSET");

class Ice
{
public:
  Ice(int argc, char** argv);

private:
  void saveWave(const CBuffer2D& wave);
  void tem();

  MultisliceParameters mParams;
};

Ice::Ice(int argc, char** argv)
{
  std::cerr << "Ice Client\n";

  if (verbosity >= 1)
    std::cout << buildString() << '\n';

  // parse the command line
  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
  {
    std::cout << "Simulation started\n";

    // output the command line options
    std::cout << CL::configString() << '\n';
  }

  if (Mode::tem == mode)
    tem();
}

void Ice::tem()
{
  if (verbosity >= 1)
    std::cout << "Inelastic mode\n";

  // load the sample description
  auto sample = std::make_shared<Sample>(sampleName);

  for (auto& atom : sample->atoms())
    atom.setBFactor(bFactor);

  // set the simulation parameters
  mParams = multisliceParamsFromCommandLine();

  InelasticSimulation inelasticSimulation(mParams, numRepetitions, sample);

  if (isTerminal(std::cout))
    inelasticSimulation.onProgress = progress;

  // execute the simulation
  inelasticSimulation.run();

  RBuffer2D temp = inelasticSimulation.mcf().diagonal();
  for (int64_t j = 0; j < rows; ++j)
    for (int64_t i = 0; i < cols; ++i)
      temp.pixel(i, j) = std::sqrt(temp.pixel(i, j));

  temp.save(outName + ".tiff", mParams, "");
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Ice ice(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
