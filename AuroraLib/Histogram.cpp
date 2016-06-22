//--- Aurora/AuroraLib/Histogram.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Histogram.hpp"

#include <fstream>

namespace Aurora
{

void Histogram::saveAsText(const std::string& filename,
                           Real trueMean, Real trueStddev) const
{
  std::ofstream outFile(filename);

  outFile << "# Mean: " << mean() << "\n"
             "# Standard deviation: " << sqrt(variance()) << "\n"
             "# Mode: " << mode() << "\n";

  if (!std::isnan(trueMean))
    outFile << "# True mean: " << trueMean << "\n";

  if (!std::isnan(trueStddev))
    outFile << "# True standard deviation: " << trueStddev << '\n';

  for (size_t i = 0; i < numBins(); ++i)
    outFile << i << ' ' << binToValue(i) << ' ' << mBins[i] << '\n';
}

} // namespace Aurora

