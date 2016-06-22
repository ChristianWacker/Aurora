//--- Aurora/Clients/Compare/Compare.cpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Physics.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <iostream>

namespace Aurora
{

enum class Mode
{
  complex, real, absSqr
};

CL::Option<Mode> mode("mode", "", Mode::complex,
  enumVal(Mode::complex), enumVal(Mode::real), enumVal(Mode::absSqr));

CL::Option<std::string> inName1 ("i1", "Directory", "NOTSET");
CL::Option<std::string> inName2 ("i2", "Directory", "NOTSET");
CL::Option<std::string> outName ("o" , "Directory", "NOTSET");

class Subtract
{
public:
  Subtract(int argc, char** argv);
};

Subtract::Subtract(int argc, char** argv)
{
  std::cerr << "Subtract Client\n";

  if (!CL::parse(argc, argv, "Subtract"))
    return;

  if (verbosity >= 1)
    std::cout << buildString() << '\n';

  std::string description;

  Buffer2DInfo info1;
  CBuffer2D buffer1(inName1, info1, description);

  Buffer2DInfo info2;
  CBuffer2D buffer2(inName2, info2, description);

  if (info1 != info2)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "Buffers are incompatible\n";
    return;
  }

  const int64_t rows = buffer1.rows();
  const int64_t cols = buffer1.cols();

  CBuffer2D bufferDiff(cols, rows);

  if (Mode::complex == mode)
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        bufferDiff.pixel(x, y) = buffer1.pixel(x, y) -
                                 buffer2.pixel(x, y);
  }
  else if (Mode::real == mode)
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        bufferDiff.pixel(x, y) = Complex(real(buffer1.pixel(x, y)) -
                                         real(buffer2.pixel(x, y)));
  }
  else if (Mode::absSqr == mode)
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        bufferDiff.pixel(x, y) = Complex(absSqr(buffer1.pixel(x, y)) -
                                         absSqr(buffer2.pixel(x, y)));
  }

  bufferDiff.save(outName, info1, description);
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Subtract subtract(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
