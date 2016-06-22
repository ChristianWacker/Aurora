//--- Aurora/AuroraLib/MultisliceOptions.hpp -----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MULTISLICE_OPTIONS_HPP
#define AURORA_AURORA_LIB_MULTISLICE_OPTIONS_HPP

#include "AuroraLib/Prerequisites.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/MultisliceParameters.hpp"

namespace Aurora
{

CL::Option<ThermicScattering> thermic("thermic", "", ThermicScattering::bFactor,
  enumVal(ThermicScattering::bFactor), enumVal(ThermicScattering::phonon));
CL::Option<std::string> sampleName
  ("sample", "Name of the sample file", "NOTSET");
CL::Option<Real> bFactor("bFactor", "", R(0.0));
CL::Option<Bandlimit> bandlimit("bandlimit", "", Bandlimit::Smooth,
  enumVal(Bandlimit::Disabled), enumVal(Bandlimit::Sharp),
  enumVal(Bandlimit::Smooth));
CL::Option<int64_t> cols("cols", "", 128);
CL::Option<int64_t> rows("rows", "", 128);
CL::Option<Real> cutoff
  ("cutoff", "Maximal extension of the atomic potentials", R(5.0));
CL::Option<Real> deltaZ
  ("deltaZ", "Step size in z direction", R(1.0));
CL::Option<Real> energy
  ("energy", "Energy of the) incomming electron in eV", R(1e5));
CL::Option<Real> minimumX("minimumX", "", R(0.0));
CL::Option<Real> minimumY("minimumY", "", R(0.0));
CL::Option<Real> minimumZ("minimumZ", "", R(0.0));
CL::Option<Real> maximumX("maximumX", "", R(1.0));
CL::Option<Real> maximumY("maximumY", "", R(1.0));
CL::Option<Real> maximumZ("maximumZ", "", R(1.0));
CL::Option<size_t> numRepetitions("numRepetitions", "", 0);
CL::Option<std::string> formFactorName("formFactor", "", "Peng");
CL::Option<std::string> phaseShiftName("phaseShift", "", "Opt");
CL::Option<std::string> potentialName("potential", "", "Opt");

inline MultisliceParameters multisliceParamsFromCommandLine()
{
  MultisliceParameters result;

  result.setThermicScattering(thermic);
  result.setResolution(cols, rows);
  result.setBox(Vector3R(minimumX, minimumY, minimumZ),
                Vector3R(maximumX, maximumY, maximumZ),
                deltaZ);
  result.setBandlimit(bandlimit);
  result.setCutoff(cutoff);
  result.setEnergy(energy);
  result.setFormFactorName(formFactorName);
  result.setPhaseShiftName(phaseShiftName);
  result.setPotentialName(potentialName);

  return result;
}

} // namespace Aurora

#endif
