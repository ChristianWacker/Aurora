//--- Aurora/AuroraLib/StemSimulation.hpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_STEM_SIMULATION_HPP
#define AURORA_AURORA_LIB_STEM_SIMULATION_HPP

#include "AuroraLib/MicroscopeSimulation.hpp"
#include "AuroraLib/StemProbe.hpp"

namespace Aurora
{

/// Simulation of a scanning transmission electron microscope (STEM)
class AURORA_API StemSimulation : public MicroscopeSimulation
{
public:
  StemSimulation(const MultisliceParameters& params,
                 size_t numRepetitions, const ConstSamplePtr& sample,
                 const std::shared_ptr<CtfCoherent>& ctf,
                 const std::string& propagatorName,
                 const std::shared_ptr<StemDetector>& detector);

  /// start the simulation
  void run() override;

  CBuffer2D& result()
  {
    return mStemResult;
  }

  typedef std::function<void (const CBuffer2D&)> OnStemResult;
  OnStemResult onStemResult;

private:
  // The wave function inside the sample
  CBuffer2D mWave;
  CBuffer2D mStemResult;
  std::string mPropagatorName;
  std::shared_ptr<StemDetector> mDetector;

  StemProbe mProbe;
};

/// Simulate the creation of an image in a scanning electron micoscope (SEM).
/// An image series is created based on a single input file. The information
/// is increased line by line. A red spot repressent the current scan
/// position.
AURORA_API void simulateStem(const std::string& inName,
                             const std::string& outName);

} // namespace Aurora

#endif
