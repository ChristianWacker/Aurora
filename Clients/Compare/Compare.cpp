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

CL::Option<std::string> inDirectoyName("i", "Directory", "NOTSET");
CL::Option<Energy> energy("energy", "Energy of the electron in eV", 80e3);
CL::Option<int64_t> index("index", "", 90);

class Compare
{
public:
  Compare(int argc, char** argv);

private:
  void compare(Real k, int64_t i);
};

Compare::Compare(int argc, char** argv)
{
  std::cerr << "Compare Client\n";

  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
    std::cerr << buildString() << '\n';

  const Real k = energy.wavenumber();

  std::cerr << "Wavelength = " << energy.wavelength() << " Angstrom\n";

  std::string params;

  Buffer2DInfo info1;
  CBuffer2D buffer1(inDirectoyName + "/1st-" + std::to_string(index) + ".tiff",
                    info1, params);

  Buffer2DInfo info2;
  CBuffer2D buffer2(inDirectoyName + "/2nd-" + std::to_string(index) + ".tiff",
                    info2, params);

  if (buffer1.cols() != buffer2.cols() || buffer2.rows() != buffer2.rows())
  {
    std::cout << "Buffer sizes do not match!\n";
    return;
  }

  CBuffer2D result(buffer1.cols(), buffer2.rows());

  Real mean = 0.0;
  Real maximum = 0.0;
  Real variance = 0.0;
  for (int64_t y = 0; y < buffer1.rows(); ++y)
  {
    for (int64_t x = 0; x < buffer1.cols(); ++x)
    {
      Real temp = 2 * k * abs(buffer1.pixel(x, y)) / abs(buffer2.pixel(x, y));
      result.pixel(x, y) = Complex(temp);
      mean += temp;
      maximum = std::max(maximum, temp);
    }
  }

  mean /= buffer1.rows() * buffer1.cols();

  for (int64_t y = 0; y < buffer1.rows(); ++y)
    for (int64_t x = 0; x < buffer1.cols(); ++x)
      variance += powerOf<2>(result.pixel(x, y).real() - mean);

  variance /= buffer1.rows() * buffer1.cols();

  std::cout << index << ' ' << mean << ' ' << sqrt(variance) << ' '
            << maximum << '\n';

  result.save("result" + std::to_string(index) + ".tiff", info1, "");
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Compare compare(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
