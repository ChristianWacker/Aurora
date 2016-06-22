//--- Aurora/Clients/OpSys/OpSys.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Ctf.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Physics.hpp"
#include "AuroraLib/Quadrature.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>

namespace Aurora
{

enum class Mode
{
  elastic, eftem, phonon, diffPhonon
};

CL::Option<Mode> mode("mode", "", Mode::elastic,
  enumVal(Mode::elastic), enumVal(Mode::eftem), enumVal(Mode::phonon),
  enumVal(Mode::diffPhonon));

CL::Option<std::string> inFilename("i", "Name of the input file");
CL::Option<std::string> outName("o", "Name of the output file");
CL::Option<std::string> ctfFilename("ctf", "Name of the output CTF file");

CL::Option<Real> cc("CC", "C_C in Angstrom", 1.5e7);
CL::Option<size_t> hermiteOrder("hermiteOrder", "", 10);
CL::Option<Real> energySpread("energySpread", "Energy spread in eV", 0.5);
CL::Option<Real> energy("energy", "Energy of the incomming electron in eV", R(1e5));

#if 0
CL::Option<Real> highTensionRipple("highTensionRipple", "", 1e-6);
CL::Option<Real> lenyCurrentInstability("lenyCurrentInstability", "", 1e-6);
#endif

class OpSys
{
public:
  OpSys(int argc, char** argv);

private:
  void elastic();
  void eftem();
  void phonon();
  void diffPhonon();
};

OpSys::OpSys(int argc, char** argv)
{
  std::cerr << "OpSys Client\n";

  if (!CL::parse(argc, argv, "Simulation of the Optical System"))
    return;

  std::cout << CL::configString() << '\n';

  if (verbosity >= 1)
    std::cout << buildString() << '\n';

  if (Mode::elastic == mode)
    elastic();
  else if (Mode::eftem == mode)
    eftem();
  else if (Mode::phonon == mode)
    phonon();
  else if (Mode::diffPhonon == mode)
    diffPhonon();
}

void OpSys::elastic()
{
  if (verbosity >= 1)
    std::cout << "Elastic mode\n";

  // load the input file
  Buffer2DInfo bufferInfo;
  std::string description;
  CBuffer2D waveIn(inFilename, bufferInfo, description);

  //std::cout << "Wavelength = " << wavelength << " Angstrom\n";

  auto waveFourier = waveIn.fft(Direction::forward);

  const int64_t cols = bufferInfo.cols();
  const int64_t rows = bufferInfo.rows();

  const Real width  = bufferInfo.width();
  const Real height = bufferInfo.height();

  const Real energy = bufferInfo.energy();

  auto ctf = CtfCoherent::fromCommandLine(energy);

  std::cout << "Kinetic energy = " << ctf->energy() * R(1e-3) << " keV\n";
  auto ctfBuffer = ctf->buffer(Buffer2DInfoBase(cols, rows, width, height),
                               R(1.0) / (rows * cols));

  if (!(ctfFilename.empty()))
  {
    // todo: reciprocal width and height
    Buffer2DInfo ctfBufferInfo(cols, rows, 0.0, 0.0, 0.0);
    ctfBuffer.save(ctfFilename, ctfBufferInfo, "Contrast transfer function");
  }

  waveFourier *= ctfBuffer;
  auto outWave = waveFourier.fft(Direction::backward);
  outWave.save(outName, bufferInfo, description);
}

void OpSys::eftem()
{
  if (verbosity >= 1)
    std::cout << "EFTEM mode\n";

  // load the metadata
  Buffer2DInfo bufferInfo;
  std::string description;
  CBuffer2D waveIn(inFilename + "0.tiff", bufferInfo, description);
  const int64_t cols = bufferInfo.cols();
  const int64_t rows = bufferInfo.rows();

  CBuffer2D result(cols, rows, Complex(0.0));

  Real focusSpread = cc * std::sqrt(powerOf<2>(energySpread / energy));

  HermitePolynomial polynomial(hermiteOrder);
  for (size_t i = 0; i < hermiteOrder; ++i)
  {
    Real focusShift = polynomial.root(i) * Math::sqrt2 * focusSpread;

    // load the input files
    CBuffer2D wave;
    wave.load(inFilename + std::to_string(i) + ".tiff",
              bufferInfo, description);

    //std::cout << "Wavelength = " << wavelength << " Angstrom\n";

    auto waveFourier = wave.fft(Direction::forward);

    const Real width  = bufferInfo.width();
    const Real height = bufferInfo.height();

    const Real bufferEnergy = bufferInfo.energy();

    std::cout << "Focus shift = " << focusShift << " Angstrom\n";

    std::cout << "Buffer energy = " << bufferEnergy << " eV\n";
    Real energyShift = focusShift / cc * energy;
    std::cout << "Energy = " << energy + energyShift << " eV\n";

    auto ctf = CtfCoherent::fromCommandLine(energy + energyShift);

    auto ctfBuffer = ctf->buffer(Buffer2DInfoBase(cols, rows, width, height),
                                 R(1.0) / (rows * cols), focusShift);

    waveFourier *= ctfBuffer;

    wave = waveFourier.fft(Direction::backward);

    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        wave.pixel(x, y) = absSqr(wave.pixel(x, y));

    result += polynomial.weight(i) / std::sqrt(Math::pi) * wave;
  }

  result.save(outName, bufferInfo, description);
}

void OpSys::phonon()
{
  if (verbosity >= 1)
    std::cout << "Brightfield phonon mode\n";

  // load the metadata
  Buffer2DInfo bufferInfo;
  std::string description;
  CBuffer2D waveIn(inFilename + "0.tiff", bufferInfo, description);

  const int64_t cols = bufferInfo.cols();
  const int64_t rows = bufferInfo.rows();
  const Real width  = bufferInfo.width();
  const Real height = bufferInfo.height();

  CBuffer2D result(cols, rows, Complex(0.0));

  auto ctf = CtfCoherent::fromCommandLine(bufferInfo.energy());
  auto ctfBuffer = ctf->buffer(Buffer2DInfoBase(cols, rows, width, height),
                               R(1.0) / (rows * cols));

  int numPhononPasses = 0;
  bool ok = true;
  do
  {
    try
    {
      // load the input files
      CBuffer2D wave(inFilename + std::to_string(numPhononPasses) + ".tiff",
                     bufferInfo, description);

      auto waveFourier = wave.fft(Direction::forward);

      waveFourier *= ctfBuffer;
      wave = waveFourier.fft(Direction::backward);

      #pragma omp parallel for
      for (int64_t y = 0; y < rows; ++y)
        for (int64_t x = 0; x < cols; ++x)
          wave.pixel(x, y) = absSqr(wave.pixel(x, y));

      result += wave;
    }
    catch (Exception&)
    {
      ok = false;
    }
  } while(ok);

  std::cout << numPhononPasses << " phonon passes.\n";

  result *= R(1.0) / numPhononPasses;
  result.save(outName, bufferInfo, description);
}

void OpSys::diffPhonon()
{
  if (verbosity >= 1)
    std::cout << "Diffraction phonon mode\n";

  // load the metadata
  Buffer2DInfo bufferInfo;
  std::string description;
  CBuffer2D waveIn(inFilename + "0.tiff", bufferInfo, description);

  const int64_t cols = bufferInfo.cols();
  const int64_t rows = bufferInfo.rows();

  // the result will be calculated in Fourier space
  CBuffer2D resultFourier(cols, rows, Complex(0.0));

  int numPhononPasses = 0;
  bool ok = true;
  do
  {
    try
    {
      // load the input files
      CBuffer2D wave(inFilename + std::to_string(numPhononPasses) + ".tiff",
                     bufferInfo, description);

      auto waveFourier = wave.fft(Direction::forward);

      const Real bufferEnergy = bufferInfo.energy();
      std::cout << "Buffer energy = " << bufferEnergy << " eV\n";

      #pragma omp parallel for
      for (int64_t y = 0; y < rows; ++y)
        for (int64_t x = 0; x < cols; ++x)
          resultFourier.pixel(x, y) += absSqr(waveFourier.pixel(x, y));

      ++numPhononPasses;
    }
    catch (Exception&)
    {
      ok = false;
    }

  } while(ok);

  std::cout << numPhononPasses << " phonon passes.\n";

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      resultFourier.pixel(x, y) = std::sqrt(resultFourier.pixel(x, y));

  resultFourier *= 1 / (bufferInfo.numPixels<Real>() * std::sqrt(static_cast<Real>(numPhononPasses)));
  resultFourier.fft(Direction::backward).save(outName, bufferInfo,
                                              description);
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::OpSys opsys(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
