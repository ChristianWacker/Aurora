//--- Aurora/Clients/Reduce/Reduce.cpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/MatrixSupport.hpp"
#include "AuroraLib/Physics.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <iostream>

namespace Aurora
{

enum class Mode
{
  abs, absSqr, current, real, scatteredCurrent
};

CL::Option<Mode> mode("mode", "", Mode::absSqr,
  enumVal(Mode::abs), enumVal(Mode::absSqr), enumVal(Mode::current),
  enumVal(Mode::real), enumVal(Mode::scatteredCurrent));
CL::Option<std::string> inName("i", "", homeDirectory());
CL::Option<std::string> inName2("i2", "", homeDirectory());
CL::Option<std::string> xCoord("xCoord", "", "0");

class Reduce
{
public:
  Reduce(int argc, char** argv);

private:
  void abs();
  void absSqr();
  void current();
  void real();
  void scatteredCurrent();
};

Reduce::Reduce(int argc, char** argv)
{
  std::cerr << "Reduce Client\n";

  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
    std::cerr << buildString() << '\n';

  if (Mode::abs == mode)
    abs();
  else if (Mode::absSqr == mode)
    absSqr();
  else if (Mode::current == mode)
    current();
  else if (Mode::real == mode)
    real();
  else if (Mode::scatteredCurrent == mode)
    scatteredCurrent();
}

void Reduce::abs()
{
  std::string params;
  Buffer2DInfo info;
  CBuffer2D wave(inName, info, params);

  std::cout << xCoord << ' ' << std::fixed << std::setprecision(12)
    << wave.absReduce() / wave.numPixels() << '\n';
}

void Reduce::absSqr()
{
  std::string params;

  Buffer2DInfo info;
  CBuffer2D wave(inName, info, params);

  std::cout << xCoord << ' ' << (wave.absSqrReduce() / wave.numPixels()) << '\n';
}

void Reduce::current()
{
  std::string params;

  Buffer2DInfo info;
  CBuffer2D wave(inName, info, params);

  VectorNC waveVector;
  bufferToVector(wave, waveVector);

  MatrixNC u;
  VectorNR eigenValuesSqr;
  loadEVDecomposition(inName2, eigenValuesSqr, u);

  VectorNC eigenValues;
  sqrtEigenValues(eigenValuesSqr, eigenValues);

  VectorNC v = u.adjoint() * waveVector;

  Real result = Real();

  #pragma omp parallel for reduction (+:result)
  for (int64_t i = 0; i < v.size(); ++i)
    result += Aurora::absSqr(v[i]) * eigenValues[i].real();

  std::cout << xCoord << ' '
            << std::setprecision(12) << result / wave.numPixels()
            << '\n';
}

void Reduce::scatteredCurrent()
{
  Buffer2DInfo info;
  std::string params;
  CBuffer2D wave(inName, info, params);

  CBuffer2D waveFft = wave.fft(Direction::forward);

  // undo the fft scaling
  #pragma omp parallel for
  for (int64_t j = 0; j < info.rows(); ++j)
  {
    // coordinates in reciprocal space
    const Real kY = j > info.rows() / 2 ? j - info.rows() : j;

    for (int64_t i = 0; i < info.cols(); ++i)
    {
      const Real kX = i > info.cols() / 2 ? i - info.cols() : i;
      // spatial frequency squared
      const Real kPerpSqr = kX * kX + kY * kY;

      if (kPerpSqr < 4)
        waveFft.pixel(i, j) = 0;
      else
        waveFft.pixel(i, j) /= info.numPixels<Real>(); // undo the fft scaling
    }
  }

  wave = waveFft.fft(Direction::backward);

  VectorNC waveVector;
  bufferToVector(wave, waveVector);

  MatrixNC u;
  VectorNR eigenValuesSqr;
  loadEVDecomposition(inName2, eigenValuesSqr, u);

  VectorNC eigenValues;
  sqrtEigenValues(eigenValuesSqr, eigenValues);

  VectorNC v = u.adjoint() * waveVector;

  Real result = 0;

  #pragma omp parallel for reduction (+:result)
  for (int64_t i = 0; i < v.size(); ++i)
    result += Aurora::absSqr(v[i]) * eigenValues[i].real();

  std::cout << xCoord << ' ' << std::fixed << std::setprecision(12)
            << result / wave.numPixels() << '\n';
}

void Reduce::real()
{
  std::string params;
  Buffer2DInfo info;
  CBuffer2D wave(inName, info, params);

  std::cout << xCoord << ' ' << std::fixed << std::setprecision(12)
    << wave.realReduce() / wave.numPixels() << '\n';
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Reduce reduce(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
