//--- Aurora/Clients/Reslice/Reslice.cpp ---------------------------------------
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

CL::Option<std::string> inName ("i", "Input Template", "NOTSET");
CL::Option<int64_t> numSlices("numSlices", "Number of Slices", 1);
CL::Option<int64_t> row("row", "Row for slicing", -1);
CL::Option<std::string> outName("o", "Output Template" , "NOTSET");

class Reslice
{
public:
  Reslice(int argc, char** argv);
};

Reslice::Reslice(int argc, char** argv)
{
  std::cerr << "Reslice Client\n";

  if (!CL::parse(argc, argv, "Slice"))
    return;

  if (verbosity >= 1)
    std::cerr << buildString() << '\n';

  Buffer2DInfo info;
  std::string description;

  std::vector<CBuffer2D> buffers(numSlices);

  buffers[0].load(inName + "0.tiff", info, description);
  for (size_t i = 1; i < static_cast<size_t>(numSlices); ++i)
  {
    Buffer2DInfo newInfo;
    buffers[i].load(inName + std::to_string(i) + ".tiff", newInfo, description);
    if (info != newInfo)
      AURORA_THROW(EInvalidParameter, "Buffers do not match");
  }

  const int64_t cols = info.cols();
  const int64_t rows = info.rows();

  if (-1 == row)
    row.setValue(rows / 2);

  std::cerr << "Row: " << row << '\n';

  if ("NOTSET" == outName)
  {
    std::cout << "# x        z        value\n";

    // output to std::out
    for (int64_t z = 0; z < numSlices; ++z)
    {
      std::cout << '\n';
      for (int64_t x = 0; x < cols; ++x)
      {
        std::cout << x << ' ' << z << ' '
                  << buffers[static_cast<size_t>(z)].pixel(x, row).real()
                  << '\n';
      }
    }
  }
  else
  {
    CBuffer2D result(cols, numSlices);
    for (int64_t z = 0; z < numSlices; ++z)
      for (int64_t x = 0; x < cols; ++x)
        result.pixel(x, z) = buffers[static_cast<size_t>(z)].pixel(x, row);

    Buffer2DInfo resultInfo(cols, numSlices, info.width(), numSlices,
                            info.energy());
    result.save(outName, resultInfo, "");
  }
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Reslice reslice(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
