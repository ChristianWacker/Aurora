//--- Aurora/Clients/Ice/InelaticSimulation.hpp --------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_ICE_INELASTIC_SIMULATION_HPP
#define AURORA_CLIENTS_ICE_INELASTIC_SIMULATION_HPP

#include "AuroraLib/MicroscopeSimulation.hpp"

#include "MutualCoherenceFunction.hpp"

namespace Aurora
{

class InelasticSimulation : public MicroscopeSimulation
{
public:
  InelasticSimulation(const MultisliceParameters& params,
                      size_t numRepetitions, const ConstSamplePtr& sample);

  void run() override;

  CBuffer2D test(Buffer2DInfo info, Real radius);

  MutualCoherenceFunction& mcf()
  {
    return mMcf;
  }

  typedef std::function<void (size_t, size_t)> OnProgress;
  OnProgress onProgress = [] (size_t, size_t) {};

private:
  void propagate();

  MutualCoherenceFunction mMcf;
  CBuffer2D mPropagatorHalf;
  CBuffer2D mPropagatorFull;
  std::unique_ptr<PhaseShift> mPhaseShift;
  size_t mNumSlices;
};

} // namespace Aurora

#endif
