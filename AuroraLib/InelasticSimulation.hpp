////////////////////////////////////////////////////////////////////////////////
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef AURORA_INELASTIC_SIMULATION_HPP
#define AURORA_INELASTIC_SIMULATION_HPP

#include "Definitions.hpp"

#include "DynArray.hpp"
#include "Math.hpp"
#include "Physics.hpp"
#include "Vector3.hpp"

#include <functional>
#include <vector>

namespace Aurora
{

class AURORA_API InelasticSimulation
{
public:
  InelasticSimulation(std::ostream& outStream);

  virtual ~InelasticSimulation() { }
  virtual void simulate(const ElasticParametersPtr& params,
                        const SamplePtr& sample);

  //typedef std::function<void (size_t, size_t,
  //                            const ConstRealBuffer2DPtr&)> OnPotential;
  std::function<void (int64_t, int64_t, const ConstComplexBuffer2DPtr&)> onWave;

  //OnPotential onPotential;

private:
  std::ostream& mOutStream;
};

} // namespace Aurora

#endif
