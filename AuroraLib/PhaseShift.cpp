//--- Aurora/AuroraLib/PhaseShift.cpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/PhaseShift.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/RBuffer2D.hpp"
#include "AuroraLib/Sample.hpp"

#include <iostream>

namespace Aurora
{

PhaseShift::PhaseShift(const MultisliceParameters& params,
                       PhaseShiftMode mode) :
  mParams(params), mMode(mode)
{}

/*static*/ PhaseShift::Creators& PhaseShift::creators()
{
  static Creators creators;
  return creators;
}

void PhaseShift::dump(const std::string& outDir) const
{
  if (verbosity >= 1)
    std::cout << "Dumping " << mParams.numSlices() << " phase shifts...";

  for (size_t i = 0; i < mParams.numSlices(); ++i)
  {
    std::string filename = outDir + "/phaseShift" + std::to_string(i) + ".tiff";
    transferFunction(i).save(filename, mParams, "");
  }

  if (verbosity >= 1)
    std::cout << "done." << std::endl;
}

const RBuffer2D& PhaseShiftRam::phaseShift(size_t index) const
{
  assert(PhaseShiftMode::phaseShift == mMode);
  assert(index < mPhaseShifts.size() && "Index out of range");
  return mPhaseShifts[index];
}

const CBuffer2D& PhaseShiftRam::transferFunction(size_t sliceIndex) const
{
  assert(PhaseShiftMode::transferFunction == mMode);
  assert(sliceIndex < mTransferFunctions.size() && "Index out of range");
  return mTransferFunctions[sliceIndex];
}

void PhaseShiftRam::setSample(const ConstSamplePtr& sample)
{
  if (!sample)
    return;

  if (verbosity >= 1)
    std::cout << "Starting phase shift calculations..." << std::endl;

  // Release the old memory
  mPhaseShifts.resize(0);
  mTransferFunctions.resize(0);

  const size_t numSlices = mParams.numSlices();
  const size_t numPixels = mParams.numPixels();
  if (verbosity >= 1)
  {
    std::cout << " Reserve memory for " << numSlices << " phaseShifts: "
              << Byte(numPixels * numSlices * sizeof(Real))
              << std::endl;
  }

  for (size_t i = 0; i < numSlices; ++i)
    mPhaseShifts.push_back(RBuffer2D(mParams.cols(), mParams.rows(), 0.0));

  using namespace std::chrono;
  auto time1 = high_resolution_clock::now();
  calculate(sample);
  auto time2 = high_resolution_clock::now();
  if (PhaseShiftMode::phaseShift == mMode)
    postProcessPhaseShifts();
  else
    postProcessTransferFunctions();
  auto time3 = high_resolution_clock::now();

  if (mParams.debug(Debug::phaseShifts))
    dump(mParams.debugDir());

#if 0
  // debug output
  // buffer for the sum of the phase shifts
  CBuffer2D phaseShift(mParams.cols(), mParams.rows());
  for (int64_t y = 0; y < mParams.rows(); ++y)
    for (int64_t x = 0; x < mParams.cols(); ++x)
      phaseShift.pixel(x, y) = 0.0;

  for (int64_t y = 0; y < mParams.rows(); ++y)
    for (int64_t x = 0; x < mParams.cols(); ++x)
      for (int64_t i = 0; i < mParams.numSlices(); ++i)
        phaseShift.pixel(x, y) += mPhaseShifts[i].pixel(x, y);

  const std::string path = homeDirectory() + "/Desktop/temp/";
  phaseShift.save(path + "phaseShift.tif", mParams, "");

  CBuffer2D transferFunc(mParams.cols(), mParams.rows());
  for (int64_t y = 0; y < mParams.rows(); ++y)
    for (int64_t x = 0; x < mParams.cols(); ++x)
      transferFunc.pixel(x, y) = fromArg(phaseShift.pixel(x, y).real());

  transferFunc.save(path + "transferFunc.tif", mParams, "");
#endif

  if (verbosity >= 1)
  {
    std::cout << "Phase shifts calculated. Time needed:\n"
                 "- Actual calculations: " << duration(time2 - time1) << "\n"
                 "- Post process: " << duration(time3 - time2) << "\n"
                 "- Total: " << duration(time3 - time1) << std::endl;
  }
}

void PhaseShiftRam::postProcessPhaseShifts()
{
  const int64_t numSlices = mParams.numSlices();

  if (Bandlimit::Disabled != mParams.bandlimit())
  {
    if (verbosity >= 1)
      std::cout << " Implement bandlimit..." << std::endl;

    RBandlimit bandlimit(mParams.bandlimit(), mParams);

    #pragma omp parallel for schedule(dynamic) firstprivate(bandlimit)
    for (int64_t i = 0; i < numSlices; ++i)
      bandlimit.apply(mPhaseShifts[i]);

    if (verbosity >= 1)
      std::cout << "finished" << std::endl;
  }
}

void PhaseShiftRam::postProcessTransferFunctions()
{
  // temporary buffer to construct the transfer functions
  CBuffer2D temp(mParams.cols(), mParams.rows());

  auto bandlimit = CBandlimit::create(mParams.bandlimit(), mParams, temp, temp);

  const size_t numSlices = mParams.numSlices();
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();

  if (verbosity >= 1)
  {
    const size_t numBytes = cols * rows * sizeof(Complex) * numSlices;
    std::cout << " Reserve memory for " << numSlices
              << " transfer functions: " << Byte(numBytes) << std::endl;
  }

  for (size_t i = 0; i < numSlices; ++i)
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        temp.pixel(x, y) = fromArg(mPhaseShifts[i].pixel(x, y));

    // delete the old buffer
    mPhaseShifts[i].free();

    if (Bandlimit::Disabled != mParams.bandlimit())
      bandlimit->apply();

    mTransferFunctions.push_back(temp.clone());
  }

  if (verbosity >= 1)
    std::cout << "finished" << std::endl;
}

PhaseShiftRef::PhaseShiftRef(const MultisliceParameters& params,
                             const ConstSamplePtr& sample,
                             PhaseShiftMode mode) :
  PhaseShiftRam(params, mode)
{
  setSample(sample);
}

void PhaseShiftRef::calculate(const ConstSamplePtr& sample)
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real cutoff = mParams.cutoff();
  const Real cutoffSqr = powerOf<2>(mParams.cutoff());
  const Vector3R minBox = mParams.minBox(0);
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  // factor for the scaling of the phase shifts
  const Real scalingFactor = mParams.lorentzFactor() / mParams.wavenumber();

  const int64_t numSlices = mParams.numSlices();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  #pragma omp parallel for schedule(dynamic)
  for (int64_t i = 0; i < numSlices; ++i)
  {
    // for all atoms
    for (auto atom : sample->atoms())
    {
      const Element& element = parametrization->element(atom.atomicNumber(),
                                                        atom.charge());
      const Vector3R atomPosition = atom.position();

      Real bFactor = 0;
      if (ThermicScattering::bFactor == mParams.thermicScattering())
        bFactor = atom.bFactor();

      const Real differenceZStart = mParams.zStart(i) - atomPosition.z();
      const Real differenceZEnd   = mParams.zEnd(i)   - atomPosition.z();

      // reject atoms that are too far away
      if ((differenceZStart > cutoff) || (-differenceZEnd > cutoff))
        continue;

      for (int64_t y = 0; y < rows; ++y)
      {
        const Real yPixel = y * deltaY + minBox.y();
        const Real differenceYSqr = powerOf<2>(yPixel - atomPosition.y());

        if (differenceYSqr > cutoffSqr)
          continue; // next line

        for (int64_t x = 0; x < cols; ++x)
        {
          const Real xPixel = x * deltaX + minBox.x();
          const Real differenceXSqr = powerOf<2>(xPixel - atomPosition.x());

          if (differenceXSqr > cutoffSqr)
            continue; // next pixel of the current line

          const Real rho = std::sqrt(differenceXSqr + differenceYSqr);

          // for the relevant slices
          mPhaseShifts[i].pixel(x, y) += scalingFactor *
            (element.phaseShift(rho, differenceZEnd, bFactor) -
             element.phaseShift(rho, differenceZStart, bFactor));
        }
      }
    }
  }
}

/// The class represents a lookup table for the phase shift calculations.
/// As the calculations of the phase shift can involve many operations that are
/// expensive in computational time, look up tables can speed up the evaluation.
/// This function uses a linear interpolation between the data points.
class PhaseShiftLookup
{
public:
  PhaseShiftLookup(const Element& element, Real scalingFactor,
                   Real bFactor, Real maxRho, Real maxZ,
                   size_t numDataPointsRho, size_t numDataPointsZ) :
    mMaxRho(maxRho), maxZ(maxZ),
    mNumDataPointsRho(numDataPointsRho),
    mNumDataPointsZ(numDataPointsZ),
    mTable(numDataPointsRho * numDataPointsZ)
  {
    #pragma omp parallel for
    for (size_t i = 0; i < mNumDataPointsRho; ++i)
    {
      const Real rho = static_cast<Real>(i) / (numDataPointsRho - 1) * mMaxRho;

      for (size_t j = 0; j < mNumDataPointsZ; ++j)
      {
        const Real z = static_cast<Real>(j) / (numDataPointsZ - 1) * maxZ;

        mTable[mNumDataPointsRho * j + i] =
          scalingFactor * element.phaseShift(rho, z, bFactor);
      }
    }
  }

  Real get(Real rho, Real z) const
  {
    assert(rho >= 0);

    Real zIntegral;
    const Real zFractional =
      std::modf(std::fabs(z) / maxZ * (mNumDataPointsZ - 1), &zIntegral);
    const uint32_t zInt = static_cast<uint32_t>(zIntegral);

    Real rhoIntegral;
    const Real rhoFractional =
      std::modf(rho / mMaxRho * (mNumDataPointsRho - 1), &rhoIntegral);
    const uint32_t rhoInt = static_cast<uint32_t>(rhoIntegral);

    assert(rhoInt < mNumDataPointsRho - 1);
    assert(zInt < mNumDataPointsZ - 1);

    // use bilinear interpolation for the final result
    const Real a = lerp(mTable[mNumDataPointsRho * zInt + rhoInt],
                        mTable[mNumDataPointsRho * (zInt + 1) + rhoInt],
                        zFractional);
    const Real b = lerp(mTable[mNumDataPointsRho * zInt + rhoInt + 1],
                        mTable[mNumDataPointsRho * (zInt + 1) + rhoInt + 1],
                        zFractional);

    return std::copysign(lerp(a, b, rhoFractional), z);
  }

private:
  const Real mMaxRho;
  const Real maxZ;
  const size_t mNumDataPointsRho;
  const size_t mNumDataPointsZ;

  std::vector<Real> mTable;
};

PhaseShiftOpt::PhaseShiftOpt(const MultisliceParameters& params,
                             const ConstSamplePtr& sample,
                             PhaseShiftMode mode) :
  PhaseShiftRam(params, mode)
{
  setSample(sample);
}

void PhaseShiftOpt::calculate(const ConstSamplePtr& sample)
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real cutoff = mParams.cutoff();
  const Vector3R minBox = mParams.minBox(0);
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();
  const size_t numSlices = mParams.numSlices();

  // factor for the scaling of the phase shifts
  const Real scalingFactor = mParams.lorentzFactor() / mParams.wavenumber();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  typedef std::pair<size_t,  Real> AtomKey;
  std::map<AtomKey, std::unique_ptr<PhaseShiftLookup>> lookupMap;

  // for all atoms
  for (auto atom : sample->atoms())
  {
    const Vector3R atomPosition = atom.position();

    // all in front?
    if (atomPosition.z() + cutoff < mParams.zStart(0))
      // no overlap
      continue;

    size_t firstSlice = 0;
    for (; firstSlice < numSlices; ++firstSlice)
    {
      if (atomPosition.z() - cutoff < mParams.zEnd(firstSlice))
        break;
    }

    // all behind?
    if (numSlices == firstSlice)
      // no overlap
      continue;

    size_t lastSlice = firstSlice;
    for (; lastSlice < numSlices - 1; ++lastSlice)
    {
      if (atomPosition.z() + cutoff < mParams.zEnd(lastSlice))
        break;
    }

    // Check, if we have already calculated a lookup table for this
    // combination of atomic number and b-factor. If not, calculate it now.
    Real bFactor = 0;
    if (ThermicScattering::bFactor == mParams.thermicScattering())
      bFactor = atom.bFactor();

    AtomKey key(atom.atomicNumber(), bFactor);
    auto it = lookupMap.find(key);
    if (lookupMap.end() == it)
    {
      const Element& element = parametrization->element(atom.atomicNumber(),
                                                        atom.charge());
      // todo: check multiplicators
      auto lookup = std::make_unique<PhaseShiftLookup>(element, scalingFactor,
        bFactor, mParams.cutoff() * R(1.5) * std::sqrt(R(2.0)),
        (mParams.cutoff() + mParams.deltaZ()) * R(1.5), 256, 128);

      it = lookupMap.insert(std::make_pair(key, std::move(lookup))).first;
    }

    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
    {
      Real yPixel = y * deltaY + minBox.y();
      Real differenceY = yPixel - atomPosition.y();

      if (std::fabs(differenceY) > cutoff)
        // next line
        continue;

      for (int64_t x = 0; x < cols; ++x)
      {
        Real xPixel = x * deltaX + minBox.x();
        Real differenceX = xPixel - atomPosition.x();

        if (std::fabs(differenceX) > cutoff)
          // next pixel
          continue;

        Real rho = std::hypot(differenceX, differenceY);

        // for the relevant slices
        for (size_t i = firstSlice; i <= lastSlice; ++i)
        {
          mPhaseShifts[i].pixel(x, y) +=
            it->second->get(rho, mParams.zEnd(i)   - atomPosition.z()) -
            it->second->get(rho, mParams.zStart(i) - atomPosition.z());
        }
      }
    }
  }
}

PhaseShiftOpt2::PhaseShiftOpt2(const MultisliceParameters& params,
                               const ConstSamplePtr& sample,
                               PhaseShiftMode mode) :
  PhaseShiftRam(params, mode)
{
  setSample(sample);
}

void PhaseShiftOpt2::calculate(const ConstSamplePtr& sample)
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real cutoff = mParams.cutoff();
  const Real cutoffSqr = powerOf<2>(mParams.cutoff());
  const Vector3R minBox = mParams.minBox();
  const Vector3R maxBox = mParams.maxBox();
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  // factor for the scaling of the phase shifts
  const Real scalingFactor = mParams.lorentzFactor() / mParams.wavenumber();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  typedef std::pair<size_t,  Real> AtomKey;
  std::map<AtomKey, std::unique_ptr<PhaseShiftLookup>> lookupMap;

  // for all atoms
  for (auto atom : sample->atoms())
  {
    // Check, if we have already calculated a lookup table for this
    // combination of atomic number and b-factor. If not, calculate it now.
    Real bFactor = 0;
    if (ThermicScattering::bFactor == mParams.thermicScattering())
      bFactor = atom.bFactor();

    AtomKey key(atom.atomicNumber(), bFactor);
    auto it = lookupMap.find(key);
    if (lookupMap.end() == it)
    {
      const Element& element = parametrization->element(atom.atomicNumber(),
                                                        atom.charge());
      // todo: check multiplicators
      auto lookup = std::make_unique<PhaseShiftLookup>(element, scalingFactor,
        bFactor, mParams.cutoff() * R(1.5) * std::sqrt(R(2.0)),
        (mParams.cutoff() + mParams.deltaZ()) * R(1.5), 256, 128);

      it = lookupMap.insert(std::make_pair(key, std::move(lookup))).first;
    }
  }

  const size_t n = mParams.numSlices();
  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < n; ++i)
  {
    // for all atoms
    for (auto atom : sample->atoms())
    {
      const Vector3R atomPosition = atom.position();

      const Real differenceZStart = mParams.zStart(i) - atomPosition.z();
      const Real differenceZEnd   = mParams.zEnd(i)   - atomPosition.z();

      // reject atoms that are too far away
      if ((differenceZStart > cutoff) || (-differenceZEnd > cutoff))
        continue; // next atom

      // fast reject
      if ((atomPosition.x() + cutoff < minBox.x()) ||
          (atomPosition.x() - cutoff > maxBox.x()))
        continue;

      Real bFactor = 0;
      if (ThermicScattering::bFactor == mParams.thermicScattering())
        bFactor = atom.bFactor();
      AtomKey key(atom.atomicNumber(), bFactor);
      const auto& lookup = lookupMap[key];

      for (int64_t y = 0; y < rows; ++y)
      {
        const Real yPixel = y * deltaY + minBox.y();
        const Real differenceYSqr = powerOf<2>(yPixel - atomPosition.y());

        if (differenceYSqr > cutoffSqr)
          continue; // next line

        for (int64_t x = 0; x < cols; ++x)
        {
          const Real xPixel = x * deltaX + minBox.x();
          const Real differenceXSqr = powerOf<2>(xPixel - atomPosition.x());

          if (differenceXSqr > cutoffSqr)
            continue; // next line

          const Real rho = std::sqrt(differenceXSqr + differenceYSqr);

          mPhaseShifts[i].pixel(x, y) +=
            lookup->get(rho, differenceZEnd) -
            lookup->get(rho, differenceZStart);
        }
      }
    }
  }
}

PhaseShiftFourier::PhaseShiftFourier(const MultisliceParameters& params,
                                     const ConstSamplePtr& sample,
                                     PhaseShiftMode mode) :
  PhaseShiftRam(params, mode)
{
  setSample(sample);
}

void PhaseShiftFourier::calculate(const ConstSamplePtr& sample)
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real cutoff = mParams.cutoff();
  const Vector3R minBox = mParams.minBox();
  const Vector3R maxBox = mParams.maxBox();
  const Real deltaKX = mParams.deltaKX();
  const Real deltaKY = mParams.deltaKY();

  const Real scalingFactor =
    deltaKX * deltaKY * mParams.lorentzFactor() / mParams.wavenumber();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  #pragma omp parallel for
  for (size_t i = 0; i < mParams.numSlices(); ++i)
  {
    CBuffer2D bufferFourier(cols, rows, 0);

    // for all atoms
    for (auto atom : sample->atoms())
    {
      const Vector3R atomPosition = atom.position();

      const Real differenceZStart = mParams.zStart(i) - atomPosition.z();
      const Real differenceZEnd   = mParams.zEnd(i)   - atomPosition.z();

      if (differenceZStart > cutoff)
        continue;

      if (-differenceZEnd > cutoff)
        continue;

      // shift the positions a little bit so they are only inside one unit cell
      const Real epsilon = R(1e-3);
      if ((atomPosition.x() + epsilon <  minBox.x()) ||
          (atomPosition.y() + epsilon <  minBox.y()) ||
          (atomPosition.x() + epsilon >= maxBox.x()) ||
          (atomPosition.y() + epsilon >= maxBox.y()))
        continue;

      const Element& element =
        parametrization->element(atom.atomicNumber(), atom.charge());

      Real bFactor = 0;
      if (ThermicScattering::bFactor == mParams.thermicScattering())
        bFactor = atom.bFactor();

      for (int64_t y = 0; y < rows; ++y)
      {
        const Real kY = deltaKY * (y > rows / 2 ? y - rows : y);
        const Real kYSqr = kY * kY;

        for (int64_t x = 0; x < cols; ++x)
        {
          const Real kX = deltaKX * (x > cols / 2 ? x - cols : x);
          const Real kXSqr = kX * kX;

          const Real phase = kX * (atomPosition.x() - minBox.x()) +
                             kY * (atomPosition.y() - minBox.y());

          bufferFourier.pixel(x, y) += std::polar(scalingFactor, phase) *
           (element.hybridPhaseShift(kXSqr + kYSqr, differenceZEnd, bFactor) -
            element.hybridPhaseShift(kXSqr + kYSqr, differenceZStart, bFactor));
        }
      }
    }

    CBuffer2D buffer(cols, rows, 0);
    buffer.assign(bufferFourier.fft(Direction::backward));
    mPhaseShifts[i].assign(bufferReal(buffer));
  }
}

PhaseShiftJitRef::PhaseShiftJitRef(const MultisliceParameters& params,
                                   const ConstSamplePtr& sample,
                                   PhaseShiftMode mode) :
  PhaseShift(params, mode), mSample(sample),
  mPhaseShift(mParams.cols(), mParams.rows()),
  mRBandlimit(mParams.bandlimit(), mParams)
{
  // FIXME: RBandlimit is always calculated
  if (PhaseShiftMode::transferFunction == mode)
  {
    mTransferFunction = CBuffer2D(mParams.cols(), mParams.rows());
    mCBandlimit = CBandlimit::create(mParams.bandlimit(), mParams,
                                     mTransferFunction, mTransferFunction);
  }
}

void PhaseShiftJitRef::setSample(const ConstSamplePtr& sample)
{
  mSample = sample;
}

const RBuffer2D& PhaseShiftJitRef::phaseShift(size_t sliceIndex) const
{
  assert(mMode == PhaseShiftMode::phaseShift);
  assert(sliceIndex < mParams.numSlices() && "Index out of range");

  calculate(sliceIndex);

  if (Bandlimit::Disabled != mParams.bandlimit())
    mRBandlimit.apply(mPhaseShift);

  return mPhaseShift;
}

const CBuffer2D& PhaseShiftJitRef::transferFunction(size_t sliceIndex) const
{
  assert(mMode == PhaseShiftMode::transferFunction);
  assert(sliceIndex < mParams.numSlices() && "Index out of range");

  calculate(sliceIndex);

  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      mTransferFunction.pixel(x, y) = fromArg(mPhaseShift.pixel(x, y));

  if (Bandlimit::Disabled != mParams.bandlimit())
    mCBandlimit->apply();

  return mTransferFunction;
}

void PhaseShiftJitRef::calculate(size_t sliceIndex) const
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real cutoff = mParams.cutoff();
  const Vector3R minBox = mParams.minBox(0);
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  // scaling factor for the integral over the potential
  const Real scalingFactor = mParams.lorentzFactor() / mParams.wavenumber();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  mPhaseShift.setZero();
  // for all atoms
  for (auto atom : mSample->atoms())
  {
    const Element& element = parametrization->element(atom.atomicNumber(),
                                                      atom.charge());
    const Vector3R atomPosition(atom.position());

    Real bFactor = 0;
    if (ThermicScattering::bFactor == mParams.thermicScattering())
      bFactor = atom.bFactor();

    const Real differenceZStart = mParams.zStart(sliceIndex) - atomPosition.z();
    const Real differenceZEnd   = mParams.zEnd(sliceIndex)   - atomPosition.z();

    if (differenceZStart > cutoff)
      continue;

    if (-differenceZEnd > cutoff)
      continue;

    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
    {
      const Real yPixel = y * deltaY + minBox.y();
      const Real differenceY = yPixel - atomPosition.y();

      if (std::fabs(differenceY) > cutoff)
        continue; // next line

      for (int64_t x = 0; x < cols; ++x)
      {
        const Real xPixel = x * deltaX + minBox.x();
        const Real differenceX = xPixel - atomPosition.x();

        if (std::fabs(differenceX) > cutoff)
          continue; // next pixel of the current line

        const Real rho = std::hypot(differenceX, differenceY);

        // for the relevant slices
        mPhaseShift.pixel(x, y) += scalingFactor *
          (element.phaseShift(rho, differenceZEnd  , bFactor) -
           element.phaseShift(rho, differenceZStart, bFactor));
      }
    }
  }
}

PhaseShiftJitOpt::PhaseShiftJitOpt(const MultisliceParameters& params,
                                   const ConstSamplePtr& sample,
                                   PhaseShiftMode mode) :
  PhaseShift(params, mode),
  mRBandlimit(mParams.bandlimit(), mParams)
{
  // FIXME: RBandlimit is always calculated
  for (size_t i = 0; i < mParams.jitNumSlices(); ++i)
  {
    // reserve memory for every phase shift
    mPhaseShifts.push_back(RBuffer2D(mParams.cols(), mParams.rows()));

    // if needed, reserve memory for every transfer function
    if (PhaseShiftMode::transferFunction == mode)
      mTransferFunctions.push_back(CBuffer2D(mParams.cols(), mParams.rows()));
  }

  if (PhaseShiftMode::transferFunction == mode)
  {
    mTemp = CBuffer2D(mParams.cols(), mParams.rows());
    mCBandlimit = CBandlimit::create(mParams.bandlimit(), mParams,
                                     mTemp, mTemp);
  }

  setSample(sample);
}
const RBuffer2D& PhaseShiftJitOpt::phaseShift(size_t sliceIndex) const
{
  assert(mMode == PhaseShiftMode::phaseShift);
  assert(sliceIndex < mParams.numSlices() && "Index out of range");

  if ((sliceIndex < mBeginSlice) || (sliceIndex >= mEndSlice))
    calculate(sliceIndex);

  return mPhaseShifts[sliceIndex - mBeginSlice];
}

const CBuffer2D& PhaseShiftJitOpt::transferFunction(size_t sliceIndex) const
{
  assert(mMode == PhaseShiftMode::transferFunction);
  assert(sliceIndex < mParams.numSlices() && "Index out of range");

  if ((sliceIndex < mBeginSlice) || (sliceIndex >= mEndSlice))
    calculate(sliceIndex);

  return mTransferFunctions[sliceIndex - mBeginSlice];
}

void PhaseShiftJitOpt::setSample(const ConstSamplePtr& sample)
{
  mSample = sample;

  // no phase shifts are calculated yet
  mBeginSlice = std::numeric_limits<size_t>::max();
  mEndSlice = std::numeric_limits<size_t>::max();

  // scaling factor for the integral over the potential
  const Real scalingFactor = mParams.lorentzFactor() / mParams.wavenumber();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  // prepare the lookup tables
  const Real pixelSize = std::min(mParams.deltaX(), mParams.deltaY());
  // todo: check multiplicators
  const Real maxRho = mParams.cutoff() * std::sqrt(R(2.0)) * R(1.1);
  const Real maxZ   = (mParams.cutoff() + mParams.deltaZ()) * R(1.1);

  mLookupMap.clear();
  for (auto atom : mSample->atoms())
  {
    // Check, if we have already calculated a lookup table for this
    // combination of atomic number and b-factor. If not, calculate it now.
    Real bFactor = 0;
    if (ThermicScattering::bFactor == mParams.thermicScattering())
      bFactor = atom.bFactor();

    AtomKey key(atom.atomicNumber(), bFactor);
    auto it = mLookupMap.find(key);
    if (mLookupMap.end() == it)
    {
      const Element& element = parametrization->element(atom.atomicNumber(),
                                                        atom.charge());

      auto lookup = std::make_unique<PhaseShiftLookup>(element, scalingFactor,
        bFactor, maxRho, maxZ, 2 * maxRho / pixelSize,
        2 * maxZ / mParams.deltaZ());

      it = mLookupMap.insert(std::make_pair(key, std::move(lookup))).first;
    }
  }
}

void PhaseShiftJitOpt::calculate(size_t beginSlice) const
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real cutoff = mParams.cutoff();
  const Real cutoffSqr = powerOf<2>(mParams.cutoff());
  const Vector3R minBox = mParams.minBox();
  const Vector3R maxBox = mParams.maxBox();
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  mBeginSlice = beginSlice;
  mEndSlice = std::min(mParams.numSlices(), mBeginSlice + mPhaseShifts.size());

  if (verbosity >= 1)
  {
    std::cout << "Calculating phase shifts for slices " << mBeginSlice
              << " to " << (mEndSlice - 1) << "... " << std::endl;
  }

  #pragma omp parallel for
  for (size_t i = mBeginSlice; i < mEndSlice; ++i)
  {
    // Clear all phaseshift
    mPhaseShifts[i - mBeginSlice].setZero();

    // for all atoms
    for (const auto& atom : mSample->atoms())
    {
      const Vector3R atomPosition = atom.position();

      const Real differenceZStart = mParams.zStart(i) - atomPosition.z();
      const Real differenceZEnd   = mParams.zEnd(i)   - atomPosition.z();

      // reject atoms that are too far away
      if ((differenceZStart > cutoff) || (-differenceZEnd > cutoff))
        continue; // next atom

      // fast reject the atom, if its x coordinate is too far outside the box
      if ((atomPosition.x() + cutoff < minBox.x()) ||
          (atomPosition.x() - cutoff > maxBox.x()))
        continue;

      Real bFactor = 0;
      if (ThermicScattering::bFactor == mParams.thermicScattering())
        bFactor = atom.bFactor();
      AtomKey key(atom.atomicNumber(), bFactor);
      const auto& lookup = mLookupMap.find(key);

      for (int64_t y = 0; y < rows; ++y)
      {
        const Real yPixel = y * deltaY + minBox.y();
        const Real differenceYSqr = powerOf<2>(yPixel - atomPosition.y());

        if (differenceYSqr > cutoffSqr)
          continue; // next line

        for (int64_t x = 0; x < cols; ++x)
        {
          const Real xPixel = x * deltaX + minBox.x();
          const Real differenceXSqr = powerOf<2>(xPixel - atomPosition.x());

          if (differenceXSqr > cutoffSqr)
            continue; // next pixel in the current line

          const Real rho = std::sqrt(differenceXSqr + differenceYSqr);

          // for the relevant slices
          mPhaseShifts[i - mBeginSlice].pixel(x, y) +=
            lookup->second->get(rho, differenceZEnd) -
            lookup->second->get(rho, differenceZStart);
        }
      }
    }

    if (PhaseShiftMode::transferFunction == mMode)
    {
      for (int64_t y = 0; y < rows; ++y)
      {
        for (int64_t x = 0; x < cols; ++x)
        {
          mTransferFunctions[i - mBeginSlice].pixel(x, y) =
            fromArg(mPhaseShifts[i - mBeginSlice].pixel(x, y));
        }
      }
    }
  }

  if (Bandlimit::Disabled != mParams.bandlimit())
  {
    if (PhaseShiftMode::phaseShift == mMode)
    {
      for (size_t i = mBeginSlice; i < mEndSlice; ++i)
        mRBandlimit.apply(mPhaseShifts[i - mBeginSlice]);
    }
    else
    {
      for (size_t i = mBeginSlice; i < mEndSlice; ++i)
      {
        mTemp.assign(mTransferFunctions[i - mBeginSlice]);
        mCBandlimit->apply();
        mTransferFunctions[i - mBeginSlice].assign(mTemp);
      }
    }
  }

  if (verbosity >= 1)
  {
    std::cout << "Phase shifts calculated. Time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

class RegisterPhaseShifts
{
public:
  RegisterPhaseShifts()
  {
    PhaseShift::registerClass<PhaseShiftRef>("Ref");
    PhaseShift::registerClass<PhaseShiftOpt>("Opt");
    PhaseShift::registerClass<PhaseShiftOpt2>("Opt2");
    PhaseShift::registerClass<PhaseShiftFourier>("Fourier");
    PhaseShift::registerClass<PhaseShiftJitRef>("JitRef");
    PhaseShift::registerClass<PhaseShiftJitOpt>("JitOpt");
  }
} registerPhaseShifts;

} // namespace Aurora
