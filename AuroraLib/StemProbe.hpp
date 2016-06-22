//--- Aurora/AuroraLib/StemProbe.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_STEM_PROBE_HPP
#define AURORA_AURORA_LIB_STEM_PROBE_HPP

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/MultisliceParameters.hpp"

namespace Aurora
{

/// This class represent a STEM probe.
/// @details All the optical properties of the probe are controlled by the CTF
///  object. including the aperture.
/// The probe is normalized to one. This will appear in AuroraView as
/// 1 / numPixels. This is ok, as the final stem image will be calculated using
/// numPixels probes.
class AURORA_API StemProbe
{
public:
  /// Creates a new object of the calcutation of STEM probes
  /// @param ctf
  ///  The CTF object controls the optical parameters of the probe.
  /// @param params
  ///  The PropagaterParameters structure contains the necessary simulations
  ///  parameters beside the optical parameters.
  /// @param targetBuffer
  ///  Reference to the target buffer in which the probe will be placed
  StemProbe(const std::shared_ptr<CtfCoherent>& ctf,
            const MultisliceParameters& params,
            CBuffer2D& targetBuffer);

  /// Calculate the probe at the specified location.
  /// The coordinates are relative to the simulation box
  void generate(Real centerX, Real centerY);

private:
  std::shared_ptr<CtfCoherent> mCtf;
  MultisliceParameters mParams;
  // Buffer in Fourier space for the construction of the STEM probe.
  CBuffer2D mFourierBuffer;
  const std::unique_ptr<Fft> mFft;
  Real mScalingFactor;
};

} // namespace Aurora

#endif
