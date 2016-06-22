//--- Aurora/AuroraLib/Potential.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_POTENTIAL_HPP
#define AURORA_AURORA_LIB_POTENTIAL_HPP

#include "AuroraLib/Factory.hpp"
#include "AuroraLib/MultisliceParameters.hpp"

namespace Aurora
{

/// Base class for all potential calculations
class AURORA_API Potential :
  public Factory<Potential, const MultisliceParameters&, const ConstSamplePtr&>
{
public:
  virtual ~Potential() {}

  size_t numPotentials() const;

  /// Returns the relativistically uncorrected potential at the position of the
  /// layer with the index @p index.
  /// @remarks
  ///  The returned potential is calculated using the Debye-Waller factor. The
  ///  returned potential is bandlimited if requested.
  /// @remarks
  ///  Before this function can be called @ref calculate() must be called.
  /// @param index
  ///  The index of the layer. This index is not identical to the slice index,
  ///  as one slice consists of three layers: one at the beginning, one in the
  ///  middle and one at the end. This is a requirement of the Runge-Kutta
  ///  algorithms.
  virtual const RBuffer2D& potential(size_t index) const = 0;

  virtual void setSample(const ConstSamplePtr& sample) = 0;

  static Creators& creators();

  void dump(const std::string& outDirectory) const;

protected:
  Potential(const MultisliceParameters& params);

  MultisliceParameters mParams;
  // the z values for which we need to calculate the potential
  std::vector<Real> mZValues;
};

/// Base class for the potential classes that the keep all the data in main
/// memory.
class AURORA_API PotentialRam : public Potential
{
protected:
  using Potential::Potential;
  PotentialRam(const PotentialRam&) = delete;
  PotentialRam& operator=(const PotentialRam& other) = delete;

  void setSample(const ConstSamplePtr& sample) override;

  std::vector<RBuffer2D> mPotentials;

public:
  const RBuffer2D& potential(size_t index) const override;

private:
  /// Implement the bandlimit
  virtual void postProcess();
  /// The actual calculation of the potentials.
  virtual void calculate(const ConstSamplePtr& sample) = 0;
};

/// Reference implementation of the potential calculations
class AURORA_API PotentialRef : public PotentialRam
{
public:
  PotentialRef(const MultisliceParameters& params,
               const ConstSamplePtr& sample);

private:
  void calculate(const ConstSamplePtr& sample) override;
};

/// Optimized implementation of the potential calculations using lookup tables.
/// The implementation is parallelized on the pixels.
class AURORA_API PotentialOpt : public PotentialRam
{
public:
  PotentialOpt(const MultisliceParameters& params,
               const ConstSamplePtr& sample);

private:
  void calculate(const ConstSamplePtr& sample) override;
};

/// Optimized implementation of the potential calculations using lookup tables.
/// The implementation is parallelized on the potentials.
class AURORA_API PotentialOpt2 : public PotentialRam
{
public:
  PotentialOpt2(const MultisliceParameters& params,
                const ConstSamplePtr& sample);

private:
  void calculate(const ConstSamplePtr& sample) override;
};

class AURORA_API PotentialFourier : public PotentialRam
{
public:
  PotentialFourier(const MultisliceParameters& params,
                   const ConstSamplePtr& sample);

private:
  void calculate(const ConstSamplePtr& sample) override;
};

class PotentialLookup;

class AURORA_API PotentialJitRef : public Potential
{
public:
  PotentialJitRef(const MultisliceParameters& params,
                  const ConstSamplePtr& sample);

  const RBuffer2D& potential(size_t index) const override;
  void setSample(const ConstSamplePtr& sample) override;

private:
  ConstSamplePtr     mSample;
  mutable RBuffer2D  mPotential;
  mutable RBandlimit mBandlimit;
};

class AURORA_API PotentialJitOpt : public Potential
{
public:
  PotentialJitOpt(const MultisliceParameters& params,
                  const ConstSamplePtr& sample);
  PotentialJitOpt(const PotentialJitOpt&) = delete;

  const RBuffer2D& potential(size_t index) const override;
  void setSample(const ConstSamplePtr& sample) override;

private:
  ConstSamplePtr mSample;

  mutable std::vector<RBuffer2D>  mPotentials;
  mutable RBandlimit mBandlimit;

  mutable size_t mBeginIndex;
  mutable size_t mEndIndex;

  typedef std::pair<size_t,  Real> AtomKey;

  #if (AURORA_PLATFORM == AURORA_PLATFORM_LINUX)
  std::map<AtomKey, std::unique_ptr<PotentialLookup>> mLookupMap;
  #elif (AURORA_PLATFORM == AURORA_PLATFORM_WINDOWS)
  std::map<AtomKey, std::shared_ptr<PotentialLookup>> mLookupMap;
  #endif

  void calculate(size_t beginIndex) const;
};

#if 0
class AURORA_API PotentialFile : public PotentialMainMemory
{
public:
  PotentialFile(const MultisliceParameters& params,
                const ConstSamplePtr& sample);

  RBuffer2D& potential(size_t index) override;

private:
  void calculate(const ConstSamplePtr&) override;
};
#endif

} // namespace Aurora

#endif
