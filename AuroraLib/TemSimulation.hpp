//--- Aurora/AuroraLib/TemSimulation.hpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_TEM_SIMULATION_HPP
#define AURORA_AURORA_LIB_TEM_SIMULATION_HPP

#include "AuroraLib/MicroscopeSimulation.hpp"

#include <functional>

namespace Aurora
{

/// Simulation of a transmission electron microscope (TEM)
class AURORA_API TemSimulation : public MicroscopeSimulation
{
public:
  TemSimulation(const MultisliceParameters& params, size_t numRepetitions,
                const ConstSamplePtr& sample,
                const std::string& propagatorName);

  /// start the simulation
  void run() override;

  /// set the input wave to a plane wave in z-direction
  void setPlaneWave();
  void setCbed(Real convergenceAngle);

  CBuffer2D& wave() { return mWave; }

  typedef std::function<void (size_t, size_t, const CBuffer2D&)> OnWave;
  OnWave onWave;

private:
  std::string mPropagatorName;

  // Internal repetitions counter.
  size_t mRepetition;

  CBuffer2D mWave;

  void fireWave(size_t index, size_t numSlices);
};

} // namespace Aurora

#endif
