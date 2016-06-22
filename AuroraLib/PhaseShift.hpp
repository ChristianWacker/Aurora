//--- Aurora/AuroraLib/PhaseShift.cpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_PHASE_SHIFT_HPP
#define AURORA_AURORA_LIB_PHASE_SHIFT_HPP

#include "AuroraLib/Bandlimit.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Factory.hpp"
#include "AuroraLib/MultisliceParameters.hpp"
#include "AuroraLib/RBuffer2D.hpp"

namespace Aurora
{

enum class PhaseShiftMode
{
  phaseShift, transferFunction
};

/// Base class for the calculation of phase shifts.
class AURORA_API PhaseShift :
  public Factory<PhaseShift, const MultisliceParameters&, const ConstSamplePtr&,
                 PhaseShiftMode>
{
public:
  virtual ~PhaseShift() {}

  /// Returns the phase shift based on
  ///  @f[
  ///    \varphi(\boldsymbol{\rho},\, z_1,\, z_2)=-\frac{m\,\gamma}{k\,\hbar^2}\
  ///    \intop_{z_1}^{z_2}\text{d}z\, V(\boldsymbol{\rho},\, z)
  ///  @f]
  /// where each @f$\boldsymbol{\rho}@f$ corresponds to a pixel in the returned
  /// buffer and the integral boundaries @f$z_i@f$ and @f$z_f@f$ are determined
  /// by the @p sliceIndex. @f$k@f$ denotes the relativistically corrected
  /// wave number.
  /// @remarks
  ///  If specified the returned buffer will be bandlimited. The returned phase
  ///  shift includes the Debye-Waller factor.
  virtual const RBuffer2D& phaseShift(size_t index) const = 0;

  /// Returns the transfer function based on
  ///  @f[
  ///    t(\boldsymbol{\rho},\, z_1,\, z_2)=\exp\left(-
  ///    \frac{\text{i}\,m\,\gamma}{k\,\hbar^2}\intop_{z_i}^{z_f}\text{d}z\, \
  ///    V(\boldsymbol{\rho},\, z)\right)
  ///  @f]
  /// where each @f$\boldsymbol{\rho}@f$ corresponds to a pixel in the returned
  /// buffer and the integral boundaries @f$z_i@f$ and @f$z_f@f$ are determined
  /// by the @p sliceIndex. @f$k@f$ denotes the relativistically corrected
  /// wave number.
  /// @remarks
  ///  If specified the returned buffer will be bandlimited. The returned
  ///  transfer function includes the Debye-Waller factor.
  virtual const CBuffer2D& transferFunction(size_t sliceIndex) const = 0;

  virtual void setSample(const ConstSamplePtr& sample) = 0;

  static Creators& creators();

  void dump(const std::string& outDirectory) const;

protected:
  PhaseShift(const MultisliceParameters& params, PhaseShiftMode mode);

  MultisliceParameters mParams;
  PhaseShiftMode mMode;
};

/// Base class for all phase shift calculations on the CPU.
class AURORA_API PhaseShiftRam : public PhaseShift
{
protected:
  using PhaseShift::PhaseShift;
  PhaseShiftRam(const PhaseShiftRam&) = delete;
  PhaseShiftRam& operator=(const PhaseShiftRam&) = delete;

  void setSample(const ConstSamplePtr& sample) override;

  std::vector<RBuffer2D> mPhaseShifts;
  std::vector<CBuffer2D> mTransferFunctions;

public:
  const RBuffer2D& phaseShift(size_t sliceIndex) const override;
  const CBuffer2D& transferFunction(size_t sliceIndex) const override;

private:
  /// The actual calculation of the phase shifts. This will be used not only by
  /// @ref calculatePhaseShifts() but also by @ref calculateTransferFunctions().
  virtual void calculate(const ConstSamplePtr& sample) = 0;

  /// Implement the bandlimit
  virtual void postProcessPhaseShifts();
  virtual void postProcessTransferFunctions();
};

/// Reference implementation for the phase shift calculations
class AURORA_API PhaseShiftRef : public PhaseShiftRam
{
public:
  PhaseShiftRef(const MultisliceParameters& params,
                const ConstSamplePtr& sample, PhaseShiftMode mode);

private:
  void calculate(const ConstSamplePtr& sample) override;
};

/// Optimized phase shift calculation using lookup tables. Outer loop iterates
/// over the atoms
class AURORA_API PhaseShiftOpt : public PhaseShiftRam
{
public:
  PhaseShiftOpt(const MultisliceParameters& params,
                const ConstSamplePtr& sample, PhaseShiftMode mode);

private:
  void calculate(const ConstSamplePtr& sample) override;
};

/// Optimized phase shift calculation using lookup tables. Outer loop iterates
/// over the slices.
class AURORA_API PhaseShiftOpt2 : public PhaseShiftRam
{
public:
  PhaseShiftOpt2(const MultisliceParameters &params,
                 const ConstSamplePtr &sample, PhaseShiftMode mode);

private:
  void calculate(const ConstSamplePtr& sample) override;
};

/// Calculate the phase shift in Fourier space.
class AURORA_API PhaseShiftFourier : public PhaseShiftRam
{
public:
  PhaseShiftFourier(const MultisliceParameters& params,
                    const ConstSamplePtr& sample, PhaseShiftMode mode);

private:
  void calculate(const ConstSamplePtr& sample) override;
};

/// Class to calculate the phase shift just-in-time (JIT)
class AURORA_API PhaseShiftJitRef : public PhaseShift
{
public:
  PhaseShiftJitRef(const MultisliceParameters& params,
                   const ConstSamplePtr& sample, PhaseShiftMode mode);

  const RBuffer2D& phaseShift(size_t sliceIndex) const override;
  const CBuffer2D& transferFunction(size_t sliceIndex) const override;

  void setSample(const ConstSamplePtr& sample) override;

private:
  ConstSamplePtr mSample;
  mutable RBuffer2D mPhaseShift;

  void calculate(size_t sliceIndex) const;

  mutable CBuffer2D mTransferFunction;

  mutable RBandlimit mRBandlimit;
  std::unique_ptr<CBandlimit> mCBandlimit;
};

class PhaseShiftLookup;

/// Calculate several phase shifts just in time using a lookup table.
/// The number of phase shifts that is calculated in parallel is controlled
/// jitNumSlices() of params.
class AURORA_API PhaseShiftJitOpt : public PhaseShift
{
public:
  PhaseShiftJitOpt(const MultisliceParameters& params,
                  const ConstSamplePtr& sample, PhaseShiftMode mode);

  const RBuffer2D& phaseShift(size_t sliceIndex) const override;
  const CBuffer2D& transferFunction(size_t sliceIndex) const override;

  void setSample(const ConstSamplePtr& sample) override;

private:
  ConstSamplePtr mSample;

  virtual void calculate(size_t beginSlice) const;

  mutable std::vector<RBuffer2D> mPhaseShifts;
  mutable std::vector<CBuffer2D> mTransferFunctions;

  mutable RBandlimit mRBandlimit;
  mutable CBuffer2D mTemp;
  std::unique_ptr<CBandlimit> mCBandlimit;

  mutable size_t mBeginSlice;
  mutable size_t mEndSlice;

  typedef std::pair<size_t,  Real> AtomKey;
  #if (AURORA_COMPILER == AURORA_COMPILER_MSVC)
  std::map<AtomKey, std::shared_ptr<PhaseShiftLookup>> mLookupMap;
  #else
  std::map<AtomKey, std::unique_ptr<PhaseShiftLookup>> mLookupMap;
  #endif
};

} // namespace Aurora

#endif
