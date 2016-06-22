//--- Aurora/Clients/Ice/InelasticSimulation.cpp -------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "InelasticSimulation.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Propagator.hpp"

#include <iostream>

namespace Aurora
{

InelasticSimulation::InelasticSimulation(const MultisliceParameters& params,
                                         size_t numRepetitions,
                                         const ConstSamplePtr& sample) :
  MicroscopeSimulation(params, numRepetitions, sample),
  mPropagatorHalf(fresnelPropagator(params, params.deltaZ() / 2,
                                    params.bandlimit())),
  mPropagatorFull(fresnelPropagator(params, params.deltaZ(),
                                    params.bandlimit())),
  mPhaseShift(PhaseShift::create(params.phaseShiftName(), params, sample,
                                 PhaseShiftMode::transferFunction)),
  mNumSlices(params.numSlices())
{}

void InelasticSimulation::propagate()
{
  mMcf.transform(Direction::forward);
  mMcf.apply(mPropagatorHalf);

  for (size_t i = 0; i < mNumSlices - 1; ++i)
  {
    mMcf.transform(Direction::backward);

    mMcf.apply(mPhaseShift->transferFunction(i));

    onProgress(i, mNumSlices);

    mMcf.transform(Direction::forward);

    mMcf.apply(mPropagatorFull);
  }

  mMcf.transform(Direction::backward);

  mMcf.apply(mPhaseShift->transferFunction(mNumSlices - 1));

  mMcf.transform(Direction::forward);

  mMcf.apply(mPropagatorHalf);

  mMcf.transform(Direction::backward);
}

void InelasticSimulation::run()
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  printStartupMessages();

  // plane incomming wave
//  CBuffer2D buffer(mParams.cols(), mParams.rows(), 1.0);

  // the initial mcf
 // mMcf = MutualCoherenceFunction::pure(buffer);

  mMcf = MutualCoherenceFunction::source(mParams, 100.0);
  std::cout << "Memory for MCF: " << Byte(mMcf.numBytes()) << '\n';

  for (size_t i = 0; i <= numRepetitions(); ++i)
    propagate();

  if (verbosity >= 1)
  {
    std::cout << "Total time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

CBuffer2D InelasticSimulation::test(Buffer2DInfo info, Real radius)
{
  const int64_t cols = info.cols();
  const int64_t rows = info.rows();
  const Real radiusSqr = radius * radius;

  CBuffer2D result(cols, rows);

  const Real centerX = info.width() / 2;
  const Real centerY = info.height() / 2;

  // Fill the buffer first
  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
  {
    const Real yPos = y * info.deltaY() - centerY;
    for (int64_t x = 0; x < cols; ++x)
    {
      const Real xPos = x * info.deltaX() - centerX;
      const Real rSqr = xPos * xPos + yPos * yPos;
      if (rSqr < radiusSqr)
        result.pixel(x, y) = Complex(1.0);
      else
        result.pixel(x, y) = Complex(0.0);
    }
  }

  return result;
}

} // namespace Aurora
