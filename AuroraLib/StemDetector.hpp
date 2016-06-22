//--- Aurora/AuroraLib/StemDetector.hpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_STEM_DETECTOR_HPP
#define AURORA_AURORA_LIB_STEM_DETECTOR_HPP

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"

namespace Aurora
{

/// Base class for a STEM detector
class AURORA_API StemDetector
{
public:
  StemDetector(Buffer2DInfo info) :
    mInfo(info)
  {}

  virtual ~StemDetector() { }

  virtual Real signal() = 0;

  /// Connects the StemDetector to an input buffer. The use of fixed buffers
  /// is necessary to improve the FFT performance.
  virtual void connect(const CBuffer2D& inputBuffer);

  virtual CBuffer2D detectorMask(Real scaling) const = 0;

protected:
  Buffer2DInfo         mInfo;
  /// the Fourier transformed input buffer
  CBuffer2D            mInputBufferFft;
  std::unique_ptr<Fft> mFft;
};

/// Annular Dark Field (ADF) Detector
/// @details
///  The detector is a ring whose inner and outer diameters are given by alpha
///  (inclusive) and beta (exlusive), respectively.
/// Examples:
///                       alpha [mrad]    beta [mrad]
///  brightfield (BF)          0              10
///  ADF                      10              50
///  high-angle (HA)-ADF      50             inf
class AURORA_API AdfDetector : public StemDetector
{
public:
  AdfDetector(Buffer2DInfo info, Real alpha, Real beta);

  /// Return the signal strength on the detector
  Real signal() override;

  CBuffer2D detectorMask(Real scaling) const override;

private:
  Real                mAlpha, mBeta;
  CBuffer2D           mDetectorMask;
};

} // namespace Aurora

#endif
