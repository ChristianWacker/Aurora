//--- Aurora/Clients/Converter/Converter.cpp -----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <iostream>

namespace Aurora
{

enum class Mode
{
  real, imag, abs, absSqr, phase
};

CL::Option<Mode> mode("mode", "", Mode::absSqr,
  enumVal(Mode::real), enumVal(Mode::imag), enumVal(Mode::abs),
  enumVal(Mode::absSqr), enumVal(Mode::phase));
CL::Option<bool> fft("fft", "", false);
CL::Option<bool> log("log", "", false);
CL::Option<std::string> inFilename("i", "Name of the input file");
CL::Option<std::string> outFilename("o", "Name of the output file");
CL::Option<Real> minValue("min", "Minimal value to map", 0);
CL::Option<Real> maxValue("max", "Maximal value to map", 1);

class Converter
{
public:
  Converter(int argc, char** argv);
};

Converter::Converter(int argc, char** argv)
{
  std::cerr << "Converter Client\n";

  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
    std::cout << buildString() << '\n';

  Buffer2DInfo info;
  std::string description;
  CBuffer2D buffer(inFilename, info, description);

  // reduce form complex values to real values
  std::function<Real (const Complex&)> extract;

  Real min = minValue;
  Real max = maxValue;

  if (Mode::real == mode)
  {
    extract = Aurora::real;
  }
  else if (Mode::imag == mode)
  {
    extract = Aurora::imag;
  }
  else if (Mode::abs == mode)
  {
    extract = Aurora::abs;
  }
  else if (Mode::absSqr == mode)
  {
    extract = Aurora::absSqr;
  }
  else if (Mode::phase == mode)
  {
    extract = Aurora::arg;
    // overwrite the user defined min and max
    min = -Math::pi;
    max = +Math::pi;
  }

  const int64_t cols = buffer.cols();
  const int64_t rows = buffer.rows();

  if (fft)
    buffer = buffer.fft(Direction::forward).fftShift();

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      buffer.pixel(x, y).real(extract(buffer.pixel(x, y)));

  if (log)
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        buffer.pixel(x, y).real(std::log(std::abs(buffer.pixel(x, y))));
  }

  auto help = [min, max] (const Complex& z)
  {
    return Gray8(saturate(R(255.0) * (z.real() - min) / (max - min)));
  };

  // FIXME: Remove static cast
  auto image = Image::fromCBuffer2D(buffer,
                                    static_cast<Image::ExtractorCGray>(help));
  image.save(outFilename);
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Converter converter(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
