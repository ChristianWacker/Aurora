//--- Aurora/Clients/Matrix/Pass.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "Pass.hpp"

#include <chrono>
#include <iostream>

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Histogram.hpp"
#include "AuroraLib/Version.hpp"

namespace Aurora
{

std::string Pass::matrixName(const std::string& prefix, size_t index) const
{
  assert(!mMatrixDir.empty() && "Filename template is empty.");

  if (std::numeric_limits<size_t>::max() == index)
    return mMatrixDir + prefix + ".mat";

  return mMatrixDir + prefix + toString(index, 4) + ".mat";
}

std::string Pass::tempName(const std::string& prefix, size_t index) const
{
  assert(!mMatrixDir.empty() && "Filename template is empty.");
  return tempDirectory() + prefix + toString(index, 4) + ".mat";
}

void Pass::dumpHistogram(const MatrixNC& matrix, const std::string& filename)
{
  const int64_t numPixexls = mParams.numPixels();
  Real minValue = std::numeric_limits<Real>::max();
  Real maxValue = 0;
  Real mean = 0;

  for (int64_t j = 0; j < numPixexls; ++j)
  {
    for (int64_t i = 0; i < numPixexls; ++i)
    {
      minValue = std::min(abs(matrix(i, j)), minValue);
      maxValue = std::max(abs(matrix(i, j)), maxValue);
      mean += abs(matrix(i, j));
    }
  }

  mean /= numPixexls;

  Histogram h(minValue * R(0.9), maxValue * R(1.1), 128);
  Real variance = 0;

  for (int64_t j = 0; j < numPixexls; ++j)
  {
    for (int64_t i = 0; i < numPixexls; ++i)
    {
      h.insert(abs(matrix(i, j)));
      variance += powerOf<2>(abs(matrix(i, j)) - mean);
    }
  }
  variance /= numPixexls;

  h.saveAsText(filename, mean, std::sqrt(variance));
}

void PassWaveOut::saveWave(const CBuffer2D& wave, const Buffer2DInfo& params,
  const std::string& infix, size_t index) const

{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  std::string filename =
    mOutTemplate + '/' + infix + '-' + std::to_string(index) + ".tiff";

  if (verbosity >= 1)
    std::cout << "Saving wave \"" << filename << "\"... " << std::flush;

  std::ostringstream description;
  description << "# Matrix, Version: " << version << "\n"
                 "# " << Time::now() << '\n'
              << CL::configString() << '\n';

  wave.save(filename, params, description.str());

  if (verbosity >= 1)
  {
    std::cout << "done. Time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

PassFinal::PassFinal(const MultisliceParameters &params,
    const std::string& matrixDir, size_t power,
    const std::string& outTemplate) :
  Pass(params, matrixDir), PassWaveOut(outTemplate), mPower(power)
{ }

Part::Part(const MultisliceParameters& params, size_t numParts, size_t part) :
  mNumParts(numParts), mPart(part)
{
  assert(mNumParts >= 1 && "Invalid number of parts.");

  // calculate the first and the last slice
  if (!(1 <= mPart && mPart <= numParts))
    AURORA_THROW(ECommandLineError, "Invalid part.");

  const size_t slicesPerJob = params.numSlices() / numParts;
  if (0 == slicesPerJob)
    AURORA_THROW(EInvalidParameter, "Too many parts.");
  size_t remainder = params.numSlices() % numParts;

  mStartSlice = 0;
  mEndSlice = slicesPerJob;
  if (remainder)
  {
    --remainder;
    ++mEndSlice;
  }

  for (size_t i = 1; i < mPart; ++i)
  {
    mStartSlice = mEndSlice;
    mEndSlice += slicesPerJob;
    if (remainder)
    {
      --remainder;
      ++mEndSlice;
    }
  }

  if (verbosity >= 1)
  {
    std::cout << "Part " << mPart << " of " << numParts << " ("
                 "startSlice: " << mStartSlice << ", "
                 "endSlice: " << mEndSlice << ")\n";
  }
}

} // namespace Aurora
