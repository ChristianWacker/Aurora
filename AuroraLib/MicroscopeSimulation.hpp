//--- Aurora/AuroraLib/MicroscopeSimulation.hpp --------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MICROSCOPE_SIMULATION_HPP
#define AURORA_AURORA_LIB_MICROSCOPE_SIMULATION_HPP

#include "AuroraLib/MultisliceParameters.hpp"

namespace Aurora
{

/// Simulation of a transmission electron microscope (TEM)
class AURORA_API MicroscopeSimulation
{
public:
  MicroscopeSimulation(const MultisliceParameters& params,
                       size_t numRepetitions, const ConstSamplePtr& sample) :
    mParams(params), mNumRepetitions(numRepetitions), mSample(sample)
  {}
  virtual ~MicroscopeSimulation() { }

  /// start the simulation
  virtual void run() = 0;

  size_t numRepetitions() const { return mNumRepetitions; }
  ConstSamplePtr sample() const { return mSample; }
  ConstSamplePtr displacedSample() const;

protected:
  MultisliceParameters mParams;

  void printStartupMessages() const;

private:
  size_t         mNumRepetitions;
  ConstSamplePtr mSample;
};

} // namespace Aurora

#endif
