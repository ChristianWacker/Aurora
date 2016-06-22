//--- Aurora/AuroraLib/Progress.hpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Progress.hpp"

#include "AuroraLib/Formatter.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>

namespace Aurora
{

void progress(int64_t part, int64_t total)
{
  using namespace std::chrono;
  static high_resolution_clock::time_point startTime;
  static high_resolution_clock::time_point lastTime;
  auto currentTime = high_resolution_clock::now();

  if (0 == part)
    startTime = currentTime;

  // write always the final number, but update only every 200 ms
  if ((total - 1 != part) &&
      (duration_cast<milliseconds>(currentTime - lastTime).count() < 200))
    return;

  high_resolution_clock::duration eta;
  if (part > 0)
  {
    auto sinceStart = currentTime - startTime;
    eta = sinceStart * total / part;
    eta -= sinceStart;
  }

  lastTime = currentTime;

  if (0 != part)
    std::cout << Terminal::clearLine;

  const Real fraction = static_cast<Real>(part) / (total - 1);
  std::cout << std::fixed << std::setprecision(1) << std::setw(5)
            << fraction * 100 << " %, ETA: " << duration(eta);

  // end with a line feed
  if (total - 1 == part)
    std::cout << '\n';

  std::cout.flush();
}

} // namespace Aurora

