//--- Aurora/AuroraLib/MultisliceParameters.hpp --------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MULTISLICE_PARAMETERS_HPP
#define AURORA_AURORA_LIB_MULTISLICE_PARAMETERS_HPP

#include "AuroraLib/Bandlimit.hpp"
#include "AuroraLib/Buffer2DInfo.hpp"

#include <vector>

namespace Aurora
{

struct Debug
{
  enum Flags
  {
    none        = 0,
    phaseShifts = 1 << 0,
    potentials  = 1 << 1,
    pass1       = 1 << 2,
    pass2       = 1 << 3,
    propagator  = 1 << 4,
    stemProbe   = 1 << 5,
    wave        = 1 << 6,
    all         = 0xffffffff
  };
};

enum class ThermicScattering
{
  bFactor, phonon
};

class AURORA_API MultisliceParameters : public Buffer2DInfo
{
public:
  MultisliceParameters() :
    mBandlimit(Bandlimit::Disabled), mCutoff(R(0.0)), mDeltaZ(Real()),
    mDepth(Real()), mMinBox(R(0.0), R(0.0), R(0.0)),
    mMaxBox(R(0.0), R(0.0), R(0.0)),
    mThermicScattering(ThermicScattering::bFactor),
    mDebugFlags(Debug::none)
  {
    setJitNumSlices();
  }

  void setBox(const Vector3R& minBox, const Vector3R& maxBox, Real deltaZ);

  void setBandlimit(Bandlimit bandlimit)
  {
    mBandlimit = bandlimit;
  }

  Bandlimit bandlimit() const
  {
    return mBandlimit;
  }

  void setCutoff(Real cutoff)
  {
    mCutoff = cutoff;
  }

  Real cutoff() const
  {
    return mCutoff;
  }

  bool debug(Debug::Flags flags) const
  {
    return 0 != (mDebugFlags & flags);
  }

  std::string debugDir() const
  {
    return mDebugDir;
  }

  Real deltaZ() const
  {
    return mDeltaZ;
  }

  Vector3R minBox(size_t sliceIndex) const
  {
    assert((sliceIndex < numSlices()));
    return mMinBoxes[sliceIndex];
  }

  size_t numSlices() const
  {
    return mMinBoxes.size() - 1;
  }

  size_t jitNumSlices() const
  {
    return mJitNumSlices;
  }

  void setDebug(Debug::Flags mode)
  {
    mDebugFlags |= mode;
  }

  void setDebugDir(const std::string& debugDir)
  {
    mDebugDir = debugDir;
  }

  std::string formFactorName() const
  {
    return mFormFactorName;
  }

  void setFormFactorName(const std::string& formFactorName)
  {
    mFormFactorName = formFactorName;
  }

  std::string phaseShiftName() const
  {
    return mPhaseShiftName;
  }

  void setPhaseShiftName(const std::string& phaseShiftName)
  {
    mPhaseShiftName = phaseShiftName;
  }

  std::string potentialName() const
  {
    return mPotentialName;
  }

  void setPotentialName(const std::string& potentialName)
  {
    mPotentialName = potentialName;
  }

  void setJitNumSlices(int64_t numSlices = -1);

  void unsetDebug(Debug::Flags mode)
  {
    mDebugFlags &= ~mode;
  }

  /// Returns the dimension of the simulation box in z direction in Angstrom
  Real depth() const
  {
    return mDepth;
  }

  /// Returns the z coordinate of the beginning of the slice with index
  /// @p sliceIndex.
  Real zStart(size_t sliceIndex) const
  {
    assert(sliceIndex < numSlices());
    return mMinBoxes[sliceIndex].z();
  }

  /// Returns the z coordinate of the end of the slice with index @p sliceIndex.
  Real zEnd(size_t sliceIndex) const
  {
    assert(sliceIndex < numSlices());
    return mMinBoxes[sliceIndex + 1].z();
  }

  Real zValues(size_t index) const
  {
    assert(index <= numSlices());
    return mMinBoxes[index].z();
  }

  Vector3R minBox() const
  {
    return mMinBox;
  }

  Vector3R maxBox() const
  {
    return mMaxBox;
  }

  ThermicScattering thermicScattering() const
  {
    return mThermicScattering;
  }

  void setThermicScattering(ThermicScattering thermicScattering)
  {
    mThermicScattering = thermicScattering;
  }

private:
  Bandlimit mBandlimit;
  Real mCutoff;
  Real mDeltaZ;
  Real mDepth;
  Vector3R mMinBox;
  Vector3R mMaxBox;
  std::string mFormFactorName;
  std::string mPhaseShiftName;
  std::string mPotentialName;
  std::vector<Vector3R> mMinBoxes;
  ThermicScattering mThermicScattering;
  uint32_t mDebugFlags;
  std::string mDebugDir;
  size_t mJitNumSlices;
};

} // namespace Aurora

#endif
