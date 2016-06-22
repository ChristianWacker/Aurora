//--- Aurora/AuroraLib/StemSimulation.cpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/StemSimulation.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Propagator.hpp"
#include "AuroraLib/StemDetector.hpp"

#include <chrono>
#include <iostream>

namespace Aurora
{

StemSimulation::StemSimulation(const MultisliceParameters& params,
    size_t numRepetitions, const ConstSamplePtr& sample,
    const std::shared_ptr<CtfCoherent>& ctf, const std::string& propagatorName,
    const std::shared_ptr<StemDetector>& detector) :
  MicroscopeSimulation(params, numRepetitions, sample),
  onStemResult([] (const CBuffer2D&) {}),
  mWave(mParams.cols(), mParams.rows()),
  mStemResult(mParams.cols(), mParams.rows()), mPropagatorName(propagatorName),
  mDetector(detector), mProbe(ctf, params, mWave)
{ }

void StemSimulation::run()
{
  using namespace std::chrono;
  auto startTime(high_resolution_clock::now());

  printStartupMessages();

  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  // In this buffer the final result will be constructed
  mStemResult.setZero();

  std::unique_ptr<Propagator> propagator;

  if (ThermicScattering::bFactor == mParams.thermicScattering())
  {
    propagator = Propagator::create(mPropagatorName, mWave, mParams, sample());
  }
  else
  {
    propagator = Propagator::create(mPropagatorName, mWave, mParams,
                                    displacedSample());
  }

  mDetector->connect(mWave);

  if (verbosity >= 1)
    std::cout << "Begin wave propagation..." << std::endl;

  // move the probe over the simulation area
  for (int64_t y = 0; y < rows; ++y)
  {
    if (verbosity >= 2)
      std::cout << "STEM-Line " << y << std::endl;

    for (int64_t x = 0; x < cols; ++x)
    {
      // update the wave with the new probe position
      mProbe.generate(x * deltaX, y * deltaY);

      if (mParams.debug(Debug::stemProbe))
      {
        mWave.save(mParams.debugDir() + "/" + std::to_string(y * cols + x)
                   + ".tiff", mParams, "");
      }

      // propagate the wave through the sample
      for (size_t i = 0; i <= numRepetitions(); ++i)
        propagator->propagate();

      // save the calculated intensity
      mStemResult.pixel(x, y) = mDetector->signal();
    }
  }

  onStemResult(mStemResult);

  if (verbosity >= 1)
  {
    std::cout << "Total time needed for STEM simulation: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void simulateStem(const std::string& inName, const std::string& outName)
{
  Image src(inName);

  const int64_t cols = src.cols();
  const int64_t rows = src.rows();

  // create an RGBA8 image
  Image dst(cols, rows, Image::PixelFormat::RGBA8);

  // Convert the input image from gray value to RGBA8.
  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      dst.pixel<Rgba8>(x, y) = Rgba8(src.pixel<uint8_t>(x, y));

  // remove pixels starting from the lower-right corner. Operate on 2x2 pixel
  // blocks.
  for (int64_t y = rows - 2; y >= 0; y -= 2)
  {
    for (int64_t x = cols - 2; x >= 0; x -= 2)
    {
      // Add a 2x2 red pixel block
      dst.pixel<Rgba8>(x    , y    ) = Rgba8(255, 0, 0, 255);
      dst.pixel<Rgba8>(x + 1, y    ) = Rgba8(255, 0, 0, 255);
      dst.pixel<Rgba8>(x    , y + 1) = Rgba8(255, 0, 0, 255);
      dst.pixel<Rgba8>(x + 1, y + 1) = Rgba8(255, 0, 0, 255);

      // calculate the position of the 2x2 pixel block to the right of the red
      // pixel block
      int64_t xNext = x + 2;
      int64_t yNext = y;
      if (xNext >= cols)
      {
        // past the end of the line => go to the next line
        xNext = 0;
        yNext += 2;
      }

      if (yNext < rows)
      {
        dst.pixel<Rgba8>(xNext    , yNext    ) = Rgba8(0);
        dst.pixel<Rgba8>(xNext + 1, yNext    ) = Rgba8(0);
        dst.pixel<Rgba8>(xNext    , yNext + 1) = Rgba8(0);
        dst.pixel<Rgba8>(xNext + 1, yNext + 1) = Rgba8(0);
      }

      // save the image
      int64_t index = (y / 2) * (cols / 2) + (x / 2);
      dst.save(outName + "/stem-" + std::to_string(index) + ".tiff");
    }
  }
}

} // namespace Aurora
