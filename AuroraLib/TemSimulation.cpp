//--- Aurora/AuroraLib/TemSimulation.cpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/TemSimulation.hpp"

#include "AuroraLib/Ctf.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Propagator.hpp"
#include "AuroraLib/StemProbe.hpp"

#include <algorithm>
#include <iostream>

namespace Aurora
{

TemSimulation::TemSimulation(const MultisliceParameters& params,
                             size_t numRepetitions,
                             const ConstSamplePtr& sample,
                             const std::string& propagatorName) :
  MicroscopeSimulation(params, numRepetitions, sample),
  onWave([] (int64_t, int64_t, const CBuffer2D&) {}),
  mPropagatorName(propagatorName), mRepetition(0),
  mWave(mParams.cols(), mParams.rows())
{
  setPlaneWave();
}

void TemSimulation::run()
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  printStartupMessages();

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

  using namespace std::placeholders;
  propagator->onWave = std::bind(&TemSimulation::fireWave, this, _1, _2);

  if (verbosity >= 1)
    std::cout << "Begin wave propagation..." << std::endl;

  for (mRepetition = 0; mRepetition <= numRepetitions(); ++mRepetition)
  {
    propagator->propagate();

    if (mParams.debug(Debug::wave))
    {
      const std::string filename =
        mParams.debugDir() + "/wave" + std::to_string(mRepetition) + ".tiff";
      mWave.save(filename, mParams, "");
    }
  }

  if (verbosity >= 1)
  {
    std::cout << "Probability = " << mWave.absSqrReduce() / mParams.numPixels<Real>()
              << '\n';

    std::cout << "Total time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void TemSimulation::setPlaneWave()
{
  // set everything to 1.0
  mWave.setValue(1.0);
}

void TemSimulation::setCbed(Real convergenceAngle)
{
  auto ctf = CtfCoherent::fromCommandLine(mParams.energy());
  ctf->setA(0, Complex(-mParams.cols() * mParams.deltaX() / 2,
                       -mParams.rows() * mParams.deltaY() / 2));
  ctf->setAperture(convergenceAngle);

  mWave.assign(ctf->buffer(mParams).fft(Direction::backward));
  const Real scaleFactor =
    std::sqrt(mParams.numPixels<Real>() / mWave.absSqrReduce());
  mWave *= scaleFactor;
}

void TemSimulation::fireWave(size_t index, size_t numSlices)
{
  // convert the index from a intra slice index to a global index
  const size_t totalNumSlices = (numRepetitions() + 1) * numSlices;
  const size_t totalIndex = mRepetition * numSlices + index;
  onWave(totalIndex, totalNumSlices, mWave);
}

} // namespace Aurora
