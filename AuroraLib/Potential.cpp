//--- Aurora/AuroraLib/Potential.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Potential.hpp"

#include "AuroraLib/Bandlimit.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/RBuffer2D.hpp"
#include "AuroraLib/Sample.hpp"

#include <iostream>

namespace Aurora
{

Potential::Potential(const MultisliceParameters& params) :
  mParams(params)
{
  // For the Runke-Kutta algorithm we need the potential at the beginning,
  // the middle and the end of a slice.
  size_t i = 0;
  for (size_t numSlices = mParams.numSlices(); i < numSlices; ++i)
  {
    mZValues.push_back(mParams.zStart(i));
    mZValues.push_back(R(0.5) * (mParams.zStart(i) + mParams.zEnd(i)));
  }
  mZValues.push_back(mParams.zEnd(i - 1));

  assert(mZValues.size() == 2 * mParams.numSlices() + 1);
}

size_t Potential::numPotentials() const
{
  return mZValues.size();
}

Potential::Creators& Potential::creators()
{
  static Creators creators;
  return creators;
}

void Potential::dump(const std::string& outDir) const
{
  if (verbosity >= 1)
    std::cout << "Dumping " << numPotentials() << " potentials...";

  for (size_t i = 0; i < numPotentials(); ++i)
  {
    std::string filename = outDir + "/potential" + std::to_string(i) + ".tiff";
    potential(i).save(filename, mParams, "");
  }

  if (verbosity >= 1)
    std::cout << "done." << std::endl;
}

void PotentialRam::setSample(const ConstSamplePtr& sample)
{
  if (!sample)
    return;

  if (verbosity >= 1)
    std::cout << "Starting potential calculations...\n";

  // Release the old memory
  mPotentials.resize(0);

  // Allocate the necessary memory
  const size_t numPixels = mParams.numPixels();
  if (verbosity >= 1)
  {
    std::cout << " Reserve memory for " << numPotentials() << " potentials: "
              << Byte(numPixels * numPotentials() * sizeof(Real))
              << std::endl;
  }

  for (size_t i = 0, e = numPotentials(); i < e; ++i)
    mPotentials.push_back(RBuffer2D(mParams.cols(), mParams.rows(), 0.0));

  using namespace std::chrono;
  auto time1 = high_resolution_clock::now();
  calculate(sample);
  auto time2 = high_resolution_clock::now();

  if (mParams.debug((Debug::potentials)))
    dump(mParams.debugDir());

  postProcess();
  auto time3 = high_resolution_clock::now();

  if (verbosity >= 1)
  {
    std::cout << "Potentials calculations. Time needed:\n"
                 "- Actual calculations: " << duration(time2 - time1) << "\n"
                 "- Post process: " << duration(time3 - time2) << "\n"
                 "- Total: " << duration(time3 - time1) << std::endl;
  }
}

const RBuffer2D& PotentialRam::potential(size_t index) const
{
  assert(index < numPotentials() && "Index out of range");
  return mPotentials[index];
}

void PotentialRam::postProcess()
{
  if (verbosity >= 1)
    std::cout << " Start post process..." << std::endl;

  if (Bandlimit::Disabled != mParams.bandlimit())
  {
    RBandlimit bandlimit(mParams.bandlimit(), mParams);

    const int64_t n = numPotentials();
    #pragma omp parallel for firstprivate(bandlimit)
    for (int64_t i = 0; i < n; ++i)
      bandlimit.apply(mPotentials[i]);
  }

  if (verbosity >= 1)
    std::cout << " End post process." << std::endl;
}

PotentialRef::PotentialRef(const MultisliceParameters& params,
                           const ConstSamplePtr& sample) :
  PotentialRam(params)
{
  setSample(sample);
}

void PotentialRef::calculate(const ConstSamplePtr& sample)
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Vector3R minBox = mParams.minBox(0);
  const Real cutoffSqr = powerOf<2>(mParams.cutoff());
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  const int64_t n = numPotentials();
  #pragma omp parallel for schedule(dynamic)
  for (int64_t i = 0; i < n; ++i)
  {
    for (const auto& atom : sample->atoms())
    {
      const Element& element = parametrization->element(atom.atomicNumber(),
                                                        atom.charge());
      const Vector3R atomPosition = atom.position();

      Real bFactor = 0;
      if (ThermicScattering::bFactor == mParams.thermicScattering())
        bFactor = atom.bFactor();

      const Real differenceZSqr = powerOf<2>(mZValues[i] - atomPosition.z());

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

          mPotentials[i].pixel(x, y) += element.potential(differenceXSqr +
                                                          differenceYSqr +
                                                          differenceZSqr,
                                                          bFactor);
        }
      }
    }
  }
}

class PotentialLookup
{
public:
  PotentialLookup(const Element& element, Real bFactor, Real maxRSqr,
                  size_t numDataPoints) :
    mMaxRSqr(maxRSqr),
    mNumDataPoints(numDataPoints),
    mTable(numDataPoints)
  {
    const int64_t n = mNumDataPoints;
    #pragma omp parallel for
    for (int64_t i = 0; i < n; ++i)
    {
      Real rSqr = static_cast<Real>(i) / (numDataPoints - 1) * mMaxRSqr;
      mTable[i] = element.potential(rSqr, bFactor);
    }
  }

  Real get(Real rSqr) const
  {
    assert(rSqr >= R(0.0));

    Real rSqrIntegral;
    const Real rSqrFractional =
      std::modf(rSqr / mMaxRSqr * (mNumDataPoints - 1), &rSqrIntegral);
    const size_t rSqrInt = static_cast<size_t>(rSqrIntegral);

    assert(rSqrInt < mNumDataPoints - 1);

    // use linear interpolation
    return lerp(mTable[rSqrInt], mTable[rSqrInt + 1], rSqrFractional);
  }

private:
  const Real mMaxRSqr;
  const size_t mNumDataPoints;

  std::vector<Real> mTable;
};

PotentialOpt::PotentialOpt(const MultisliceParameters& params,
                           const ConstSamplePtr& sample) :
  PotentialRam(params)
{
  setSample(sample);
}

void PotentialOpt::calculate(const ConstSamplePtr& sample)
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Real cutoff = mParams.cutoff();
  const Real cutoffSqr = powerOf<2>(mParams.cutoff());
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();
  const Vector3R minBox = mParams.minBox(0);

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  typedef std::pair<size_t, Real> AtomKey;
  std::map<AtomKey, std::unique_ptr<PotentialLookup>> lookupMap;

  for (const auto& atom : sample->atoms())
  {
    const Vector3R& atomPosition = atom.position();

    // all in front
    if (atomPosition.z() + cutoff < mZValues[0])
      continue;

    // find the first layer that the current atom can influence
    size_t firstLayer = 0;
    for (; firstLayer < numPotentials(); ++firstLayer)
    {
      if (atomPosition.z() - cutoff < mZValues[firstLayer])
        break;
    }

    // check, if the atom can influence anything
    if (numPotentials() == firstLayer)
      continue;

    // find the last layer that the current atom can influence
    size_t lastLayer = firstLayer;
    for (; lastLayer < numPotentials() - 1; ++lastLayer)
    {
      if (atomPosition.z() + cutoff < mZValues[lastLayer])
        break;
    }

    Real bFactor = 0;
    if (ThermicScattering::bFactor == mParams.thermicScattering())
      bFactor = atom.bFactor();

    // Check, if we have already calculated a lookup table for the atomic
    // number and b-factor combination. If not, calculate it now.
    AtomKey key(atom.atomicNumber(), bFactor);
    auto it = lookupMap.find(key);

    if (lookupMap.end() == it)
    {
      const Element& element = parametrization->element(atom.atomicNumber(),
                                                        atom.charge());
      // todo: check multiplicators
      auto lookup = std::make_unique<PotentialLookup>(element, bFactor,
                      cutoffSqr * R(1.1) * R(3.0), 16394);
      it = lookupMap.insert(std::make_pair(key, std::move(lookup))).first;
    }

    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
    {
      Real yPixel = y * deltaY + minBox.y();
      Real differenceYSqr = powerOf<2>(yPixel - atomPosition.y());

      if (differenceYSqr > cutoffSqr)
        continue; // next line

      for (int64_t x = 0; x < cols; ++x)
      {
        Real xPixel = x * deltaX + minBox.x();
        Real differenceXSqr = powerOf<2>(xPixel - atomPosition.x());

        if (differenceXSqr > cutoffSqr)
            continue; // next pixel of the current line

        // for the relevant layers
        for (size_t i = firstLayer; i <= lastLayer; ++i)
        {
          Real differenceZSqr = powerOf<2>(atomPosition.z() - mZValues[i]);
          mPotentials[i].pixel(x, y) += it->second->get(differenceXSqr +
                                                        differenceYSqr +
                                                        differenceZSqr);
        }
      }
    }
  }
}

PotentialOpt2::PotentialOpt2(const MultisliceParameters& params,
                             const ConstSamplePtr& sample) :
  PotentialRam(params)
{
  setSample(sample);
}

void PotentialOpt2::calculate(const ConstSamplePtr& sample)
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Vector3R minBox = mParams.minBox();
  const Vector3R maxBox = mParams.maxBox();
  const Real cutoff = mParams.cutoff();
  const Real cutoffSqr = powerOf<2>(mParams.cutoff());
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  typedef std::pair<size_t, Real> AtomKey;
  std::map<AtomKey, std::unique_ptr<PotentialLookup>> lookupMap;

  for (const Atom& atom : sample->atoms())
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
      auto lookup = std::make_unique<PotentialLookup>(element, bFactor,
                      cutoffSqr * R(1.1) * R(3.0), 16394);
      lookupMap.insert(std::make_pair(key, std::move(lookup)));
    }
  }

  const int64_t n = numPotentials();
  #pragma omp parallel for schedule(dynamic)
  for (int64_t i = 0; i < n; ++i)
  {
    for (const Atom& atom : sample->atoms())
    {
      const Vector3R atomPosition = atom.position();

      const Real differenceZSqr = powerOf<2>(mZValues[i] - atomPosition.z());

      // reject atoms that are too far away
      if (differenceZSqr > cutoffSqr)
        continue;

      // fast reject the atom, if its x coordinate is too far outside the box
      if ((atomPosition.x() + cutoff < minBox.x()) ||
          (atomPosition.x() - cutoff > maxBox.x()))
        continue;

      // reject atoms that are too far away
      Real bFactor = 0;
      if (ThermicScattering::bFactor == mParams.thermicScattering())
        bFactor = atom.bFactor();
      AtomKey key(atom.atomicNumber(), bFactor);
      const auto& lookup = lookupMap.find(key);

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

          mPotentials[i].pixel(x, y) += lookup->second->get(differenceXSqr +
                                                            differenceYSqr +
                                                            differenceZSqr);
        }
      }
    }
  }
}

PotentialFourier::PotentialFourier(const MultisliceParameters& params,
                                   const ConstSamplePtr& sample) :
  PotentialRam(params)
{
  setSample(sample);
}

void PotentialFourier::calculate(const ConstSamplePtr& sample)
{
  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Vector3R minBox = mParams.minBox();
  const Vector3R maxBox = mParams.maxBox();
  const Real cutoff = mParams.cutoff();
  const Real deltaKX = mParams.deltaKX();
  const Real deltaKY = mParams.deltaKY();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  const Real scalingFactor = deltaKX * deltaKY * powerOf<2>(Phys::hBar) /
                               (4 * Math::piSqr * Phys::mE);

  const int64_t n = numPotentials();
  #pragma omp parallel for
  for (int64_t i = 0; i < n; ++i)
  {
    CBuffer2D bufferFourier(cols, rows, 0);

    for (const auto& atom : sample->atoms())
    {
      const Vector3R atomPosition = atom.position();

      const Real differenceZ = mZValues[i] - atomPosition.z();

      if (abs(differenceZ) > cutoff)
        continue;

      // shift the position a little bit so they are only inside one unit cell
      const Real epsilon = R(1e-3);
      if ((atomPosition.x() + epsilon <  minBox.x()) ||
          (atomPosition.y() + epsilon <  minBox.y()) ||
          (atomPosition.x() + epsilon >= maxBox.x()) ||
          (atomPosition.y() + epsilon >= maxBox.y()))
        continue;

      const Element& element = parametrization->element(atom.atomicNumber(),
                                                        atom.charge());

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

          bufferFourier.pixel(x, y) -= std::polar(scalingFactor, phase) *
           element.hybridPotential(kXSqr + kYSqr, differenceZ, bFactor);
        }
      }
    }

    CBuffer2D buffer(cols, rows);
    buffer.assign(bufferFourier.fft(Direction::backward));
    mPotentials[i].assign(bufferReal(buffer));
  }
}

PotentialJitRef::PotentialJitRef(const MultisliceParameters& params,
                           const ConstSamplePtr& sample) :
  Potential(params), mSample(sample),
  mPotential(mParams.cols(), mParams.rows()),
  mBandlimit(mParams.bandlimit(), mParams)
{}

void PotentialJitRef::setSample(const ConstSamplePtr& sample)
{
  mSample = sample;
}

const RBuffer2D& PotentialJitRef::potential(size_t index) const
{
  assert(index < numPotentials() && "Index out of range");

  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Vector3R minBox = mParams.minBox(0);
  const Real cutoffSqr = powerOf<2>(mParams.cutoff());
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  mPotential.setZero();
  for (const auto& atom : mSample->atoms())
  {
    const Element& element = parametrization->element(atom.atomicNumber(),
                                                      atom.charge());
    const Vector3R atomPosition = atom.position();

    const Real differenceZSqr = powerOf<2>(mZValues[index] - atomPosition.z());

    #pragma omp parallel for
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

        mPotential.pixel(x, y) += element.potential(differenceXSqr +
                                                    differenceYSqr +
                                                    differenceZSqr);
      }
    }
  }

  if (Bandlimit::Disabled != mParams.bandlimit())
    mBandlimit.apply(mPotential);

  return mPotential;
}

PotentialJitOpt::PotentialJitOpt(const MultisliceParameters& params,
                                 const ConstSamplePtr& sample) :
  Potential(params),
  mBandlimit(mParams.bandlimit(), mParams)
{
  // reserve memory for every potential
  for (size_t i = 0; i < mParams.jitNumSlices(); ++i)
    mPotentials.push_back(RBuffer2D(mParams.cols(), mParams.rows()));

  setSample(sample);
}

const RBuffer2D& PotentialJitOpt::potential(size_t index) const
{
  assert(index < numPotentials() && "Index out of range");

  if ((index < mBeginIndex) || (index >= mEndIndex))
    calculate(index);

  return mPotentials[index - mBeginIndex];
}

void PotentialJitOpt::setSample(const ConstSamplePtr &sample)
{
  mSample = sample;

  // no potentials are calculated yet
  mBeginIndex = std::numeric_limits<size_t>::max();
  mEndIndex = std::numeric_limits<size_t>::max();

  auto parametrization =
    FormFactorParametrization::create(mParams.formFactorName());

  const Real cutoffSqr = powerOf<2>(mParams.cutoff());

  // prepare the lookup tables
  const Real pixelSize = std::min(std::min(mParams.deltaX(), mParams.deltaY()),
                                  mParams.deltaZ());

  // todo: check multiplicators
  // 3 from the diagonal, 1.1 security margin
  const Real maxRSqr = cutoffSqr * R(3.0) * R(1.1);

  mLookupMap.clear();
  for (const Atom& atom : sample->atoms())
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

      auto lookup = std::make_unique<PotentialLookup>(element, bFactor,
        maxRSqr, 2 * maxRSqr / pixelSize);
      mLookupMap.insert(std::make_pair(key, std::move(lookup)));
    }
  }
}

void PotentialJitOpt::calculate(size_t beginIndex) const
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const Vector3R minBox = mParams.minBox();
  const Vector3R maxBox = mParams.maxBox();
  const Real cutoff = mParams.cutoff();
  const Real cutoffSqr = powerOf<2>(mParams.cutoff());
  const Real deltaX = mParams.deltaX();
  const Real deltaY = mParams.deltaY();

  mBeginIndex = beginIndex;
  mEndIndex = std::min(mZValues.size(), mBeginIndex + mPotentials.size());

  if (verbosity >= 1)
  {
    std::cout << "Calculating potentials for layers " << mBeginIndex
              << " to " << (mEndIndex - 1) << "... " << std::endl;
  }

  #pragma omp parallel for
  for (size_t i = mBeginIndex; i < mEndIndex; ++i)
  {
    // clear all potentials
    mPotentials[i - mBeginIndex].setZero();

    for (const Atom& atom : mSample->atoms())
    {
      const Vector3R atomPosition = atom.position();

      const Real differenceZSqr = powerOf<2>(mZValues[i] - atomPosition.z());

      // reject atoms that are too far away
      if (differenceZSqr > cutoffSqr)
        continue;

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

          mPotentials[i - mBeginIndex].pixel(x, y) +=
            lookup->second->get(differenceXSqr +
                                differenceYSqr +
                                differenceZSqr);
        }
      }
    }
  }

  if (Bandlimit::Disabled != mParams.bandlimit())
    for (size_t i = mBeginIndex; i < mEndIndex; ++i)
      mBandlimit.apply(mPotentials[i - mBeginIndex]);

  if (verbosity >= 1)
  {
    std::cout << "Potentials calculated. Time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

#if 0
void PotentialFile::calculate(const ConstSamplePtr& sample)
{
  std::ifstream inFile(homeDirectory() + "/Desktop/hydrogen/tp.xyz");

  int i = 0;
  Real minValue = std::numeric_limits<Real>::min();
  Real maxValue = std::numeric_limits<Real>::max();
  Real sum = 0.0;
  int indexX = 0;
  int indexY = 0;
  int indexZ = 0;
  int voxels = 0;

  while (!inFile.eof())
  {
    ++i;

    std::string line;
    std::getline(inFile, line);

    // skip empty lines
    if ("" == trim(line))
      continue;

    // skip comments
    if ('#' == line[0])
      continue;

    std::istringstream lineStream(line);

    Vector3R coords;
    lineStream >> coords.x();
    lineStream >> coords.y();
    lineStream >> coords.z();

    Real value;
    lineStream >> value;

    maxValue = std::max(value, maxValue);
    minValue = std::min(value, minValue);

    if ((indexX < mParams.cols()) &&
        (indexY < mParams.rows()) &&
        (indexZ < mParams.numSlices()))
    {
      sum += value;
      ++voxels;

      // todo: scaling factor?
      //value += 0.282;
      value *= -1.0;
      value *= 27.0;

      mPotentials[indexZ].pixel(indexX, indexY) = value;
    }

    ++indexX;

    if (indexX == mParams.cols() + 1)
    {
       ++indexY;
      indexX = 0;
    }

    if (indexY == mParams.rows() + 1)
    {
      ++indexZ;
      indexY = 0;
    }
  }

  std::cout << "I: " << i << '\n';

  std::cout << "Minimal Value: " << minValue << '\n';
  std::cout << "Mean Value: "    << sum / voxels << '\n';
  std::cout << "Maximal Value: " << maxValue << '\n';
}
#endif

class RegisterPotentials
{
public:
  RegisterPotentials()
  {
    Potential::registerClass<PotentialRef>("Ref");
    Potential::registerClass<PotentialOpt>("Opt");
    Potential::registerClass<PotentialOpt2>("Opt2");
    Potential::registerClass<PotentialFourier>("Fourier");
    Potential::registerClass<PotentialJitRef>("JitRef");
    Potential::registerClass<PotentialJitOpt>("JitOpt");
  }
} registerPotentials;

} // namespace Aurora
