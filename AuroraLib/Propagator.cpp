//--- Aurora/AuroraLib/Propagator.cpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Propagator.hpp"

#include "AuroraLib/Bandlimit.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Math.hpp"
#include "AuroraLib/Potential.hpp"
#include "AuroraLib/RBuffer2D.hpp"
#include "AuroraLib/Utils.hpp"

#include <iostream>

namespace Aurora
{

using namespace std::complex_literals;

CBuffer2D fresnelPropagator(const Buffer2DInfo& info, Real deltaZ,
                            Bandlimit bandlimit)
{
  const int64_t cols = info.cols();
  const int64_t rows = info.rows();

  CBuffer2D result(cols, rows);

  const Real deltaKX = info.deltaKX();
  const Real deltaKY = info.deltaKY();

  const Real prefactor = -deltaZ / (2 * info.wavenumber());

  #pragma omp parallel for
  for (int64_t j = 0; j < rows; ++j)
  {
    // coordinates in reciprocal space
    const Real kY = deltaKY * (j > rows / 2 ? j - rows : j);

    for (int64_t i = 0; i < cols; ++i)
    {
      const Real kX = deltaKX * (i > cols / 2 ? i - cols : i);
      // spatial frequency squared
      const Real kPerpSqr = kX * kX + kY * kY;
      result.pixel(i, j) = fromArg(prefactor * kPerpSqr);
    }
  }

  if (Bandlimit::Disabled == bandlimit)
  {
    // Correct the FFT scaling
    const Real invFftScaling = 1 / info.numPixels<Real>();
    result *= invFftScaling;
  }
  else
  {
    // implement the bandlimit. We do not need to explicitly incorporate the
    // FFT scaling factor, as it is already included in the masks.
    RBuffer2D mask;
    if (Bandlimit::Sharp == bandlimit)
      mask = CBandlimitSharp::mask(info);
    else if (Bandlimit::Smooth == bandlimit)
      mask = CBandlimitSmooth::mask(info);
    else
      AURORA_UNREACHABLE;

    result *= mask;
  }

  return result;
}

Propagator::Propagator(CBuffer2D& wave, const MultisliceParameters& params) :
  onWave([] (int64_t, int64_t, const CBuffer2D&) {}),
  mParams(params), mWave(wave), mZ(params.zStart(0))
{}

Propagator::Creators& Propagator::creators()
{
  static Creators creators;
  return creators;
}

void Propagator::reportDivergence(const std::string& additionalText)
{
  const int64_t limit = 3;

  ++mNumDivergences;
  if (verbosity >= 2)
  {
    if (mNumDivergences <= limit)
      std::cout << additionalText << "Series does not seem to converge.\n";
    if (mNumDivergences == limit)
      std::cout << "Further output will be suppressed.\n";
  }
}

PropagatorPhaseShift::PropagatorPhaseShift(CBuffer2D& wave,
                                           const MultisliceParameters& params,
                                           const ConstSamplePtr& sample,
                                           PhaseShiftMode mode) :
  Propagator(wave, params),
  mPhaseShift(PhaseShift::create(params.phaseShiftName(), params, sample, mode))
{}

PropagatorPotential::PropagatorPotential(CBuffer2D& wave,
                                         const MultisliceParameters& params,
                                         const ConstSamplePtr& sample) :
  Propagator(wave, params),
  mPotential(Potential::create(params.phaseShiftName(), params, sample))
{}

PropagatorClassicalFT::PropagatorClassicalFT(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPhaseShift(wave, params, sample, PhaseShiftMode::transferFunction),
  mWaveFft(wave.createCompatibleBuffer()),
  mFft(Fft::create(mWave, mWaveFft, Direction::forward, false)),
  mFftInv(Fft::create(mWaveFft, mWave, Direction::backward, false)),
  mPropagatorHalf(fresnelPropagator(mParams, params.deltaZ() / 2,
                                    mParams.bandlimit())),
  mPropagatorFull(fresnelPropagator(mParams, params.deltaZ(),
                                    mParams.bandlimit()))
{}

void PropagatorClassicalFT::propagate()
{
  // go to the Fourier space for the first time
  mFft->transform();

  // apply this Fresnel propagator
  mWaveFft *= mPropagatorHalf;

  const size_t numSlices = mParams.numSlices();
  // for all slices but the last
  for (size_t i = 0; i < numSlices - 1; ++i)
  {
    // to direct space
    mFftInv->transform();

    // apply the potential
    mWave *= mPhaseShift->transferFunction(i);

    onWave(i, numSlices, mWave);

    if (verbosity >= 2)
    {
      std::cout << "Slice " << i << ": Probability = "
                << mWave.absSqrReduce() / mParams.numPixels<Real>()
                << std::endl;
    }

    // to Fourier space
    mFft->transform();

    // apply Fresnel propagator
    mWaveFft *= mPropagatorFull;
  }

  // to direct space
  mFftInv->transform();

  // apply potential
  mWave *= mPhaseShift->transferFunction(numSlices - 1);

  // to Fourier space
  mFft->transform();

  // the slicing of Kirkland is different
  //waveFourier->multiplyAssign(propagatorFull);

  // apply the Fresnel propagator for half a slice
  mWaveFft *= mPropagatorHalf;

  // to direct space
  mFftInv->transform();

  if (verbosity >= 2)
  {
    std::cout << "Slice " << numSlices - 1 << ": Probability = "
              << mWave.absSqrReduce() / mParams.numPixels<Real>() << std::endl;
  }

  onWave(numSlices - 1, numSlices, mWave);
}

template<class Method>
PropagatorClassical<Method>::PropagatorClassical(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPhaseShift(wave, params, sample, PhaseShiftMode::transferFunction),
  mTemp1(wave.createCompatibleBuffer()),
  mTemp2(CMarginBuffer2D<Method::margin>(mParams.cols(), mParams.rows())),
  mBandlimit(CBandlimit::create(params.bandlimit(), params, mTemp1, mWave))
{}

template<class Method>
void PropagatorClassical<Method>::propagate()
{
  // This function is written to reflect the structure of
  // PropagatorClassicalFT::propagate()
  const size_t numSlices = mParams.numSlices();

  propagateOneSlice(mWave, R(0.5) * mParams.deltaZ());

  // for all slices but the last
  for (size_t i = 0; i < numSlices - 1; ++i)
  {
    // apply potential
    mWave *= mPhaseShift->transferFunction(i);

    onWave(i, numSlices, mWave);

    if (verbosity >= 2)
    {
      std::cout << "Slice " << i << ": Probability = "
                << mWave.absSqrReduce() / mParams.numPixels<Real>()
                << std::endl;
    }

    propagateOneSlice(mWave, mParams.deltaZ());
  }

  // apply potential
  mWave *= mPhaseShift->transferFunction(numSlices - 1);

  propagateOneSlice(mWave, R(0.5) * mParams.deltaZ());

  onWave(numSlices - 1, numSlices, mWave);
}

template<class Method>
void PropagatorClassical<Method>::propagateOneSlice(CBuffer2D& wave,
                                                    Real deltaZ)
{
  const Real prefactor = R(0.5) * deltaZ /
    (mParams.wavenumber() * mParams.deltaX() * mParams.deltaY());

  // we construct the new wave in mTempBuffer1. This the "one" from the
  // expansion of the exponential function
  mTemp1.assign(wave);

  int64_t iteration = 1;
  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence("Helper: ");
      break;
    }

    // copy wave to mTemp2
    mTemp2.assign(wave);

    if (wave.info().maxAbs < R(1e-7))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    mTemp2.replicateMargin();

    wave.assign(prefactor / static_cast<Real>(iteration) *
                R(1.0i) * Laplace<Method>(mTemp2));
    mTemp1 += wave;

    ++iteration;
  }

  // implement the bandlimit
  mBandlimit->apply();
}

template<class Method>
PropagatorSingle<Method>::PropagatorSingle(CBuffer2D& wave,
                                           const MultisliceParameters& params,
                                           const ConstSamplePtr& sample) :
  PropagatorPhaseShift(wave, params, sample, PhaseShiftMode::phaseShift),
  mWaveAfter(wave.createCompatibleBuffer()),
  mTemp(CMarginBuffer2D<Method::margin>(mParams.cols(), mParams.rows())),
  mBandlimit(CBandlimit::create(params.bandlimit(), params, mWaveAfter, mWave))
{}

template<class Method>
void PropagatorSingle<Method>::propagate()
{
  const size_t numSlices = mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    propagateOneSlice(mWave, mWaveAfter, mPhaseShift->phaseShift(i));

    if (verbosity >= 2)
    {
      std::cout << "Slice " << i << ": Probability = "
             << mWaveAfter.absSqrReduce() / mParams.numPixels<Real>()
             << std::endl;
    }

    // implement the band-limit
    mBandlimit->apply();

    onWave(i, numSlices, mWave);
  }
}

template<class Method>
void PropagatorSingle<Method>::propagateOneSlice(CBuffer2D& waveIn,
                                                 CBuffer2D& waveOut,
                                                 const RBuffer2D& phaseShift)
{
  const Real prefactor = R(0.5) * mParams.deltaZ() /
    (mParams.wavenumber() * mParams.deltaX() * mParams.deltaY());

  waveOut.assign(waveIn);

  int64_t iteration = 1;
  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence();
      break;
    }

    mTemp.assign(waveIn);

    if (waveIn.info().maxAbs < R(1e-7))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    mTemp.replicateMargin();

    waveIn.assign(R(1.0i) / static_cast<Real>(iteration) *
                  (phaseShift * mTemp + prefactor * Laplace<Method>(mTemp)));
    waveOut += waveIn;

    ++iteration;
  }
}

PropagatorSingleFT::PropagatorSingleFT(CBuffer2D& wave,
                                       const MultisliceParameters& params,
                                       const ConstSamplePtr& sample) :
  PropagatorPhaseShift(wave, params, sample, PhaseShiftMode::phaseShift),
  mWaveAfter(wave.createCompatibleBuffer()),
  mTemp(params.cols(), params.rows()),
  mLaplaceFT(params, mTemp, mTemp,
             Complex(R(-0.5) * params.deltaZ() / params.wavenumber())),
  mBandlimit(CBandlimit::create(params.bandlimit(), params, mWaveAfter, mWave))
{}

void PropagatorSingleFT::propagate()
{
  const size_t numSlices = mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    propagateOneSlice(mWave, mWaveAfter, mPhaseShift->phaseShift(i));

    if (verbosity >= 2)
    {
      std::cout << "Slice " << i << ": Probability = "
             << mWaveAfter.absSqrReduce() / mParams.numPixels<Real>()
             << std::endl;
    }

    // implement the band-limit
    mBandlimit->apply();

    onWave(i, numSlices, mWave);
  }
}

void PropagatorSingleFT::propagateOneSlice(CBuffer2D& waveIn,
                                           CBuffer2D& waveOut,
                                           const RBuffer2D& phaseShift)
{
  waveOut.assign(waveIn);

  int64_t iteration = 1;
  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence();
      break;
    }

    if (waveIn.info().maxAbs < R(1e-7))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    mTemp.assign(waveIn);
    mLaplaceFT.apply();

    waveIn.assign(R(1.0i) / static_cast<Real>(iteration) *
                  (phaseShift * waveIn + mTemp));
    waveOut += waveIn;

    ++iteration;
  }
}

template<class Method>
PropagatorPartial<Method>::PropagatorPartial(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPhaseShift(wave, params, sample, PhaseShiftMode::phaseShift),
  mWaveAfter(wave.createCompatibleBuffer()),
  mHelp(CMarginBuffer2D<Method::margin>(mParams.cols(), mParams.rows())),
  mOperatorPower(wave.createCompatibleBuffer()),
  mBandlimit(CBandlimit::create(params.bandlimit(), params, mWaveAfter, mWave))
{}

template<class Method>
void PropagatorPartial<Method>::propagate()
{
  const size_t numSlices = mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    propagateOneSlice(mWave, mWaveAfter, mPhaseShift->phaseShift(i));

    if (verbosity >= 2)
    {
      std::cout << "Slice " << i << ": Probability = "
             << mWaveAfter.absSqrReduce() / mParams.numPixels<Real>()
             << std::endl;
    }

    // implement the band-limit
    mBandlimit->apply();

    onWave(i, numSlices, mWave);
  }
}

template<class Method>
void PropagatorPartial<Method>::propagateOneSlice(CBuffer2D& waveIn,
                                                  CBuffer2D& waveOut,
                                                  const RBuffer2D& phaseShift)
{
  waveOut.assign(waveIn);

  size_t iteration = 1;
  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence();
      break;
    }

    if (waveIn.info().maxAbs < R(1e-7))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    helper(waveIn, phaseShift);

    waveIn *= R(1.0) / static_cast<Real>(iteration);
    waveOut += waveIn;

    ++iteration;
  }
}

template<class Method>
void PropagatorPartial<Method>::helper(CBuffer2D& wave,
                                       const RBuffer2D& phaseShift)
{
  const Real wavenum = mParams.wavenumber();
  const Real laplaceFactor =
    R(1.0) / (powerOf<2>(wavenum) * mParams.deltaX() * mParams.deltaY());

  mOperatorPower.assign(wave);
  wave *= R(1.0i) * phaseShift;

  int64_t iteration = 1;

  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence("Helper: ");
      break;
    }

    // Copy the result so far into the help buffer.
    mHelp.assign(mOperatorPower);

    const Real binom = binomial(R(0.5), iteration);
    if (mOperatorPower.info().maxAbs * std::fabs(binom) < R(1e-5))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    mHelp.replicateMargin();

    mOperatorPower.assign(laplaceFactor * Laplace<Method>(mHelp));
    wave += binom * mParams.deltaZ() * wavenum * R(1.0i) * mOperatorPower;

    ++iteration;
  }
}

template<class Method>
PropagatorFull<Method>::PropagatorFull(CBuffer2D& wave,
                                       const MultisliceParameters& params,
                                       const ConstSamplePtr& sample) :
  PropagatorPhaseShift(wave, params, sample, PhaseShiftMode::phaseShift),
  mWaveAfter(wave.createCompatibleBuffer()),
  mHelp(CMarginBuffer2D<Method::margin>(mParams.cols(), mParams.rows())),
  mOperatorPower(wave.createCompatibleBuffer()),
  mBandlimit(CBandlimit::create(params.bandlimit(), params, mWaveAfter, mWave))
{}

template<class Method>
void PropagatorFull<Method>::propagate()
{
  const size_t numSlices = mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    propagateOneSlice(mWave, mWaveAfter, mPhaseShift->phaseShift(i));

    if (verbosity >= 2)
    {
      std::cout << "Slice " << i << ": Probability = "
                << mWaveAfter.absSqrReduce() / mParams.numPixels<Real>()
                << std::endl;
    }

    // implement the band-limit
    mBandlimit->apply();

    onWave(i, numSlices, mWave);
  }
}

template<class Method>
void PropagatorFull<Method>::propagateOneSlice(CBuffer2D& waveIn,
                                               CBuffer2D& waveOut,
                                               const RBuffer2D& phaseShift)
{
  waveOut.assign(waveIn);

  int64_t iteration = 1;
  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence();
      break;
    }

    if (waveIn.info().maxAbs < R(1e-7))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    // calculate the argument of the exponential function
    helper(waveIn, phaseShift);

    // the faculty of the exponential function
    waveIn *= R(1.0) / static_cast<Real>(iteration);
    waveOut += waveIn;

    ++iteration;
  }
}

template<class Method>
void PropagatorFull<Method>::helper(CBuffer2D& wave,
                                    const RBuffer2D& phaseShift)
{
  const Real wavenum = mParams.wavenumber();
  const Real laplaceFactor =
    R(1.0) / (powerOf<2>(wavenum) * mParams.deltaX() * mParams.deltaY());
  const Real phaseShiftFactor = R(2.0) / (wavenum * mParams.deltaZ());

  mOperatorPower.assign(wave);
  int64_t iteration = 1;

  wave.setZero();

  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence("Helper: ");
      break;
    }

    // Copy the result so far into the help buffer.
    mHelp.assign(mOperatorPower);

    const Real binom = binomial(R(0.5), iteration);
    if (mOperatorPower.info().maxAbs * std::fabs(binom) < R(1e-5))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    mHelp.replicateMargin();

    mOperatorPower.assign(laplaceFactor * Laplace<Method>(mHelp) +
                          phaseShiftFactor * phaseShift * mHelp);
    wave += wavenum * mParams.deltaZ() * binom * R(1.0i) * mOperatorPower;

    ++iteration;
  }
}

PropagatorFullFT::PropagatorFullFT(CBuffer2D& wave,
                                   const MultisliceParameters& params,
                                   const ConstSamplePtr& sample) :
  PropagatorPhaseShift(wave, params, sample, PhaseShiftMode::phaseShift),
  mWaveAfter(wave.createCompatibleBuffer()),
  mHelp(params.cols(), params.rows()),
  mTemp(params.cols(), params.rows()),
  mLaplaceFT(params, mHelp, mTemp,
             Complex(-1 / powerOf<2>(params.wavenumber()))),
  mOperatorPower(params.cols(), params.rows()),
  mBandlimit(CBandlimit::create(params.bandlimit(), params, mWaveAfter, mWave))
{}

void PropagatorFullFT::propagate()
{
  const size_t numSlices = mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    propagateOneSlice(mWave, mWaveAfter, mPhaseShift->phaseShift(i));

    if (verbosity >= 2)
    {
      std::cout << "Slice " << i << ": Probability = "
                << mWaveAfter.absSqrReduce() / mParams.numPixels<Real>()
                << std::endl;
    }

    // implement the band-limit
    mBandlimit->apply();

    onWave(i, numSlices, mWave);
  }
}

void PropagatorFullFT::propagateOneSlice(CBuffer2D& waveIn, CBuffer2D& waveOut,
                                         const RBuffer2D& phaseShift)
{
  waveOut.assign(waveIn);

  int64_t iteration = 1;
  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence();
      break;
    }

    if (waveIn.info().maxAbs < R(1e-7))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    // calculate the argument of the exponential function
    helper(waveIn, phaseShift);

    // the factorial of the exponential function
    waveIn *= R(1.0) / static_cast<Real>(iteration);
    waveOut += waveIn;

    ++iteration;
  }
}

void PropagatorFullFT::helper(CBuffer2D& wave, const RBuffer2D& phaseShift)
{
  const Real wavenum = mParams.wavenumber();
  const Real phaseShiftFactor = R(2.0) / (wavenum * mParams.deltaZ());

  mOperatorPower.assign(wave);
  int64_t iteration = 1;

  wave.setZero();

  while (true)
  {
    if (40 == iteration)
    {
      reportDivergence("Helper: ");
      break;
    }

    // Copy the result so far into the help buffer.
    mHelp.assign(mOperatorPower);

    const Real binom = binomial(0.5, iteration);

    if (mOperatorPower.info().maxAbs * std::fabs(binom) < R(1e-5))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << iteration << " iterations.\n";
      break;
    }

    mLaplaceFT.apply();

    mOperatorPower.assign(mTemp + phaseShiftFactor * phaseShift * mHelp);
    wave += wavenum * mParams.deltaZ() * binom * R(1.0i) * mOperatorPower;

    ++iteration;
  }
}

PropagatorFullFT2::PropagatorFullFT2(CBuffer2D& wave,
                                    const MultisliceParameters& params,
                                    const ConstSamplePtr& sample) :
  PropagatorPhaseShift(wave, params, sample, PhaseShiftMode::phaseShift),
  mWaveAfter(wave.createCompatibleBuffer()),
  mOperatorPower(params.cols(), params.rows()),
  mLaplaceFT(params, mWave, mOperatorPower,
             Complex(1 / powerOf<2>(params.wavenumber()))),
  mBandlimit(CBandlimit::create(params.bandlimit(), params, mWaveAfter, mWave))
{}

void PropagatorFullFT2::propagate()
{
  const size_t numSlices = mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    propagateOneSlice(mPhaseShift->phaseShift(i));

    if (verbosity >= 2)
    {
      std::cout << "Slice " << i << ": Probability = "
                << mWaveAfter.absSqrReduce() / mParams.numPixels<Real>()
                << std::endl;
    }

    // implement the band-limit
    mBandlimit->apply();

    onWave(i, numSlices, mWave);
  }
}

void PropagatorFullFT2::propagateOneSlice(const RBuffer2D& phaseShift)
{
  const Real phaseShiftFactor = -2 / (mParams.wavenumber() * mParams.deltaZ());

  mWaveAfter.assign(mWave);

  unsigned n = 1;
  while (true)
  {
    if (20 == n)
    {
      reportDivergence();
      break;
    }

    mLaplaceFT.apply();
    mOperatorPower += phaseShiftFactor * phaseShift * mWave;

    const Complex coeff = coefficient(n);

    if (mOperatorPower.info().maxAbs * std::abs(coeff) < R(1e-7))
    {
      if (verbosity >= 2)
        std::cout << "Stop after " << n << " iterations.\n";
      break;
    }

    mWaveAfter += coeff * mOperatorPower;
    mWave.assign(mOperatorPower);

    ++n;
  }
}

Complex PropagatorFullFT2::coefficient(unsigned l)
{
  Complex result = 0;

  for (unsigned m = 0; m < l; ++m)
  {
    result +=
      std::pow(R(1.0i), static_cast<Real>(l - m)) *
      std::pow(mParams.wavenumber() * mParams.deltaZ(), static_cast<Real>(l - m)) *
      Complex(factorial(l + m - 1)) /
      (std::pow(R(2.0), static_cast<Real>(m + l)) *
       factorial(l) * factorial(m) * factorial(l - m - 1));
  }

  return result;
}

template<class Method>
PropagatorRungeKutta1<Method>::PropagatorRungeKutta1(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPotential(wave, params, sample),
  mA(mParams.cols(), mParams.rows()),
  mATemp(mParams.cols(), mParams.rows()),
  mAP0(mParams.cols(), mParams.rows()),
  mAPA(mParams.cols(), mParams.rows()),
  mAPB(mParams.cols(), mParams.rows()),
  mAPC(mParams.cols(), mParams.rows()),
  mBandlimitA(CBandlimit::create(mParams.bandlimit(), mParams, wave, wave)),
  mPde(mParams)
{}

template<class Method>
void PropagatorRungeKutta1<Method>::propagate()
{
  const size_t numSlices = mParams.numSlices();

  // copy the incoming wave function in a buffer with a margin
  mA.assign(mWave);

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    const Real deltaZ = this->mParams.zEnd(i) - this->mParams.zStart(i);

    // the actual classical Runge-Kutta algorithm
    mPde.f(mA, this->mPotential->potential(2 * i), mAP0);
    mATemp.assign(mA + R(0.5) * deltaZ * mAP0);
    mPde.f(mATemp, this->mPotential->potential(2 * i + 1), mAPA);
    mATemp.assign(mA + R(0.5) * deltaZ * mAPA);
    mPde.f(mATemp, this->mPotential->potential(2 * i + 1), mAPB);
    mATemp.assign(mA + deltaZ * mAPB);
    mPde.f(mATemp, this->mPotential->potential(2 * i + 2), mAPC);

    // Combine all waves
    mA += deltaZ / R(6.0) * (mAP0 + R(2.0) * (mAPA + mAPB) + mAPC);

    if (Bandlimit::Disabled != mParams.bandlimit())
    {
      // implement the band-limit
      // copy mA to mWave, as this buffer has no margin and hence can be Fourier
      // transformed
      mWave.assign(mA);
      mBandlimitA->apply();
      mA.assign(mWave);
    }

    mWave.assign(mA);

    if (verbosity >= 2)
    {
      mPde.diagnostic(i, mZ, this->mA);
      mZ += deltaZ;
    }

    if (mParams.debug(Debug::propagator))
    {
      static int index = 0;
      if ((index < 200) || (index % 200 == 0))
      {
        Buffer2DInfo info(mParams.cols(), mParams.rows(),
                          mParams.width(), mParams.energy(), mParams.energy());
        mWave.save(mParams.debugDir() + "/A" + std::to_string(index) + ".tif",
                   info, "");
      }
      ++index;
    }

    onWave(i, numSlices, mWave);
  }
}

PropagatorRungeKutta1FT::PropagatorRungeKutta1FT(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPotential(wave, params, sample),
  mA(mParams.cols(), mParams.rows()),
  mATemp(mParams.cols(), mParams.rows()),
  mAP0(mParams.cols(), mParams.rows()),
  mAPA(mParams.cols(), mParams.rows()),
  mAPB(mParams.cols(), mParams.rows()),
  mAPC(mParams.cols(), mParams.rows()),
  mTemp(mParams.cols(), mParams.rows()),
  mLaplaceFT(mParams, mTemp, mTemp, -R(0.5i) / mParams.wavenumber()),
  mPotentialFactor(-R(1.0i) * Phys::mE * mParams.lorentzFactor() /
                   (mParams.wavenumber() * powerOf<2>(Phys::hBar))),
  mBandlimitA(CBandlimit::create(mParams.bandlimit(), mParams, mA, mA))
{}

void PropagatorRungeKutta1FT::f(CBuffer2D& inA, const RBuffer2D& potential,
                                CBuffer2D& out)
{
  mTemp.assign(inA);
  mLaplaceFT.apply();
  out.assign(mTemp + mPotentialFactor * potential * inA);
}

void PropagatorRungeKutta1FT::diagnostic(size_t slice, Real z,
                                         const CBuffer2D& a) const
{
  std::cout << "Slice " << slice << ": " << "z = " << z
            << "; Probability = "
            << a.absSqrReduce() / (mParams.cols() * mParams.rows()) << '\n';
}

void PropagatorRungeKutta1FT::propagate()
{
  const size_t numSlices = mParams.numSlices();

  // copy the incoming wave function in a buffer with a margin
  mA.assign(mWave);

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    const Real deltaZ = this->mParams.zEnd(i) - this->mParams.zStart(i);

    // the actual classical Runge-Kutta algorithm
    f(mA, mPotential->potential(2 * i), mAP0);
    mATemp.assign(mA + R(0.5) * deltaZ * mAP0);
    f(mATemp, mPotential->potential(2 * i + 1), mAPA);
    mATemp.assign(mA + R(0.5) * deltaZ * mAPA);
    f(mATemp, mPotential->potential(2 * i + 1), mAPB);
    mATemp.assign(mA + deltaZ * mAPB);
    f(mATemp, mPotential->potential(2 * i + 2), mAPC);

    // Combine all waves
    mA += deltaZ / R(6.0) * (mAP0 + R(2.0) * (mAPA + mAPB) + mAPC);

    if (Bandlimit::Disabled != mParams.bandlimit())
      mBandlimitA->apply();

    mWave.assign(mA);

    if (verbosity >= 2)
    {
      diagnostic(i, mZ, mA);
      mZ += deltaZ;
    }

    onWave(i, numSlices, mWave);
  }
}

PdeBase::PdeBase(const MultisliceParameters& params) :
  mCols(params.cols()), mRows(params.rows()), mNumPixels(mCols * mRows)
{}

inline int64_t PdeBase::index(int64_t x, int64_t y, int64_t block) const
{
  return block * mNumPixels + y * mCols + x;
}

inline int64_t PdeBase::indexS(int64_t x, int64_t y, int64_t block) const
{
  if (x < 0)
    x += mCols;
  else if (x >= mCols)
    x -= mCols;

  if (y < 0)
    y += mRows;
  else if (y >= mRows)
    y -= mRows;

  return block * mNumPixels + y * mCols + x;
}

template<class Method>
Pde1<Method>::Pde1(const MultisliceParameters& params) :
  PdeBase(params),
  mLaplaceFactor(R(0.5i) / (params.deltaX() * params.deltaY() * params.wavenumber())),
  mPotentialFactor(-R(1.0i) * Phys::mE * params.lorentzFactor() /
                   (params.wavenumber() * powerOf<2>(Phys::hBar)))
{}

template<class Method>
inline void Pde1<Method>::f(MarginBuffer& inA, const RBuffer2D& potential,
                            CBuffer2D& out) const
{
  inA.replicateMargin();
  out.assign(mLaplaceFactor * Laplace<Method>(inA) +
             mPotentialFactor * potential * inA);
}

template<class Method>
inline void Pde1<Method>::g(MarginBuffer& inA,const RBuffer2D& potential,
                            Real weight, CBuffer2D& out) const
{
  inA.replicateMargin();
  out.assign(inA + mLaplaceFactor * weight * R(1.0i) * Laplace<Method>(inA)
                 + mPotentialFactor * weight * R(1.0i) * potential * inA);
}

template<class Method>
void Pde1<Method>::diagnostic(size_t slice, Real z,
                              const MarginBuffer& a) const
{
  std::cout << "Slice " << slice << ": " << "z = " << z
            << "; Probability = " << a.absSqrReduce() / mNumPixels << std::endl;
}

template<class Method>
void Pde1<Method>::fillTriplets(Triplets& triplets, int64_t rowBlock,
                                int64_t colBlock, Real weight,
                                const Aurora::RBuffer2D& potential) const
{
  // for a detailed description, see Pde2::fill().
  const Complex weightLaplace = weight * mLaplaceFactor;
  const Complex weightPotential = weight * mPotentialFactor;

  // kronecker delta: The ones on the main diagonal from the Lobatto algorithm
  const Real kroneckerDelta = (rowBlock == colBlock) ? 1 : 0;

  for (int64_t y = 0; y < mRows; ++y)
  {
    for (int64_t x = 0; x < mCols; ++x)
    {
      const int64_t rowIndex = index(x, y, rowBlock);

      triplets.emplace_back(rowIndex, index(x, y, colBlock),
                              kroneckerDelta
                            + weightLaplace   * Laplace<Method>::factor[0]
                            + weightPotential * potential.pixel(x, y));

      // Add the non-central elements of the Laplace operator
      for (int64_t i = 1; i <= Method::margin; ++i)
      {
        const Complex factor = weightLaplace * Laplace<Method>::factor[i];
        triplets.emplace_back(rowIndex, indexS(x + i, y    , colBlock), factor);
        triplets.emplace_back(rowIndex, indexS(x - i, y    , colBlock), factor);
        triplets.emplace_back(rowIndex, indexS(x    , y + i, colBlock), factor);
        triplets.emplace_back(rowIndex, indexS(x    , y - i, colBlock), factor);
      }
    }
  }
}

template<class Method>
Pde2<Method>::Pde2(const MultisliceParameters& params) :
  PdeBase(params),
  mWaveNum(params.wavenumber()),
  mLaplaceFactor(-R(1.0) / (params.deltaX() * params.deltaY())),
  mPotentialFactor(R(2.0) * Phys::mE * params.lorentzFactor() /
                   powerOf<2>(Phys::hBar))
{}

#if 1
template<class Method>
inline void Pde2<Method>::initB(const MarginBuffer&, CBuffer2D& b,
                                const std::unique_ptr<Potential>&)
{
  b.setZero();
}
#else
template<class Method>
void Pde2<Method>::initB(const MarginBuffer& a, CBuffer2D& b,
                         const std::unique_ptr<Potential>& potential)
{
  // Calculate the first step with RungeKutta1 to calculate the corrected
  // start value for b
  const Real deltaZ = mParams.zEnd(0) - mParams.zStart(0);

  mA.assign(mWave);
  MarginBuffer mATemp(mParams.cols(), mParams.rows());
  MarginBuffer mAP0(mParams.cols(), mParams.rows());
  MarginBuffer mAPA(mParams.cols(), mParams.rows());
  MarginBuffer mAPB(mParams.cols(), mParams.rows());
  MarginBuffer mAPC(mParams.cols(), mParams.rows());

  fStart(mA, mPotential->potential(0), mAP0);
  mATemp.assign(mA + _R(0.5) * deltaZ * mAP0);
  fStart(mATemp, mPotential->potential(1), mAPA);
  mATemp.assign(mA + _R(0.5) * deltaZ * mAPA);
  fStart(mATemp, mPotential->potential(1), mAPB);
  mATemp.assign(mA + deltaZ * mAPB);
  fStart(mATemp, mPotential->potential(2), mAPC);

  // Combine all waves
  b.assign(_R(1.0) / _R(6.0) * (mAP0 + _R(2.0) * (mAPA + mAPB) + mAPC));
}
#endif

template<class Method>
inline void Pde2<Method>::f(MarginBuffer& inA, const CBuffer2D& inB,
                            const RBuffer2D& potential, CBuffer2D& outBP) const
{
  inA.replicateMargin();
  outBP.assign(- mWaveNum * R(2.0i) * inB
               + mLaplaceFactor * Laplace<Method>(inA)
               + mPotentialFactor * potential * inA);
}

template<class Method>
inline void Pde2<Method>::g(MarginBuffer& inA, const CBuffer2D& inB,
                            const RBuffer2D& potential, Real weight,
                            CBuffer2D& out) const
{
  inA.replicateMargin();
  out.assign(inB - weight * mWaveNum * R(2.0i) * inB
                 + weight * mLaplaceFactor * Laplace<Method>(inA)
                 + weight * mPotentialFactor * potential * inA);
}

template<class Method>
void Pde2<Method>::diagnostic(size_t slice, Real z, const MarginBuffer& a,
                              const CBuffer2D& b) const
{
  std::cout << "Slice " << slice << ": " << "z = " << z
         << "; mA.absSqr() = " << a.absSqrReduce() / mNumPixels
         << "; mB.absSqr() = " << b.absSqrReduce() / mNumPixels << std::endl;
}

template<class Method>
void Pde2<Method>::fillTriplets(Triplets& triplets, int64_t rowBlock,
                                int64_t columnBlock, Real weight,
                                const RBuffer2D& potential) const
{
  // We create a 2x2 block matrix that can itself be a block of another
  // blockmatrix. The position of this is block in the outer matrix is
  // determined by rowBlock and columnBlock.

  // Premultiply the weight with the scaling factors.
  const Real weightLaplace = weight * mLaplaceFactor;
  const Real weightPotential = weight * mPotentialFactor;

  // kronecker delta: This identity matrix is only present, if the current block
  // is on the main diagonal of the outer block matrix. The identity matrix is
  // required by the implicit algorithms.
  const Real kroneckerDelta = (rowBlock == columnBlock) ? 1 : 0;
  // The factor needed in the lower right block (11).
  const Complex factor11(kroneckerDelta, R(-2.0) * mWaveNum * weight);

  // We have 2 * n blocks in on dimension with numPixels rows and cols, where
  // n is the number of blocks in the outer matrix (3 for Lobatto methods).
  // Here, we calculate the absolute block indices. Remember: We construct a 2x2
  // block matrix,
  const int64_t rB0 = 2 * rowBlock + 0;
  const int64_t rB1 = 2 * rowBlock + 1;
  const int64_t cB0 = 2 * columnBlock + 0;
  const int64_t cB1 = 2 * columnBlock + 1;

  for (int64_t y = 0; y < mRows; ++y)
  {
    for (int64_t x = 0; x < mCols; ++x)
    {
      const int64_t rowIndex0 = index(x, y, rB0);
      const int64_t rowIndex1 = index(x, y, rB1);

      // block00: This is the identity matrix, iff we are on the main diagonal.
      triplets.emplace_back(rowIndex0, index(x, y, cB0), kroneckerDelta);
      // block01: This is the feedback term: B -> A
      triplets.emplace_back(rowIndex0, index(x, y, cB1), weight);

      // block10: The actual PDE.
      // This is the central term
      triplets.emplace_back(rowIndex1, index(x, y, cB0),
                              weightLaplace   * Laplace<Method>::factor[0]
                            + weightPotential * potential.pixel(x, y));

      // Add the non-central elements of the Laplace operator
      for (int64_t i = 1; i <= Method::margin; ++i)
      {
        const Real factor = weightLaplace * Laplace<Method>::factor[i];
        triplets.emplace_back(rowIndex1, indexS(x + i, y    , cB0), factor);
        triplets.emplace_back(rowIndex1, indexS(x - i, y    , cB0), factor);
        triplets.emplace_back(rowIndex1, indexS(x    , y + i, cB0), factor);
        triplets.emplace_back(rowIndex1, indexS(x    , y - i, cB0), factor);
      }

      // block 11: This is the feedback term B -> B and the other part of the
      // identity matrix.
      triplets.emplace_back(rowIndex1, index(x, y, cB1), factor11);
    }
  }
}

#if 0
template<class Method>
void Pde2<Method>::fStart(MarginBuffer& in, const RBuffer2D& potential,
                          MarginBuffer& out)
{
  const Complex factor1 = 0.5_i /
    (mParams.wavenumber() * mParams.deltaX() * mParams.deltaY());
  const Complex factor2 = -1.0_i * Phys::mE * mParams.lorentzFactor() /
    (mParams.wavenumber() * powerOf<2>(Phys::hBar));

  in.replicateMargin();
  out.assign(factor1 * Laplace<Method>(in) + factor2 * potential * in);
}
#endif

template<class Method>
Pde3<Method>::Pde3(const MultisliceParameters& params) :
  PdeBase(params),
  mWaveNum(params.wavenumber()),
  mEnergyFactor(-powerOf<2>(mWaveNum)),
  mLaplaceFactor(-R(1.0) / (params.deltaX() * params.deltaY())),
  mPotentialFactor(R(2.0) * Phys::mE * params.lorentzFactor() /
                   powerOf<2>(Phys::hBar))
{ }

template<class Method>
inline void Pde3<Method>::initB(const MarginBuffer&, CBuffer2D& b,
                                const std::unique_ptr<Potential>&)
{
  b.setValue(Complex(0.0, mWaveNum));
  mInitialBAbsSqr = b.absSqrReduce();
  std::cout << "Initial b.absSqr() = " << mInitialBAbsSqr << std::endl;
}

template<class Method>
inline void Pde3<Method>::f(MarginBuffer& inA, const CBuffer2D&,
                            const RBuffer2D& potential, CBuffer2D& outBP) const
{
  inA.replicateMargin();
  outBP.assign(  mEnergyFactor * inA
               + mLaplaceFactor * Laplace<Method>(inA)
               + mPotentialFactor * potential * inA);
}

template<class Method>
inline void Pde3<Method>::g(MarginBuffer& inA, const CBuffer2D& inB,
                            const RBuffer2D& potential, Real weight,
                            CBuffer2D& out) const
{
  inA.replicateMargin();
  out.assign(inB + weight * mEnergyFactor * inA
                 + weight * mLaplaceFactor * Laplace<Method>(inA)
                 + weight * mPotentialFactor * potential * inA);
}

template<class Method>
void Pde3<Method>::diagnostic(size_t slice, Real z, const MarginBuffer& a,
                              const CBuffer2D& b) const
{
  std::cout << "Slice " << slice << ": " << "z = " << z
            << "; mA.absSqr() = " << a.absSqrReduce() / mNumPixels
            << "; mB.absSqr() = " << b.absSqrReduce() / mInitialBAbsSqr
            << std::endl;
}

template<class Method>
void Pde3<Method>::fillTriplets(Triplets& triplets, int64_t rowBlock,
                                int64_t columnBlock, Real weight,
                                const RBuffer2D& potential) const
{
  // for a detailed description, see Pde2::fill().

  const Real weightEnergy    = weight * mEnergyFactor;
  const Real weightLaplace   = weight * mLaplaceFactor;
  const Real weightPotential = weight * mPotentialFactor;

  // kronecker delta: The ones on the main diagonal from the Lobatto algorithm
  const Real kroneckerDelta = (rowBlock == columnBlock ? 1.0 : 0.0);

  const int64_t rB0 = 2 * rowBlock + 0;
  const int64_t rB1 = 2 * rowBlock + 1;
  const int64_t cB0 = 2 * columnBlock + 0;
  const int64_t cB1 = 2 * columnBlock + 1;

  for (int64_t y = 0; y < mRows; ++y)
  {
    for (int64_t x = 0; x < mCols; ++x)
    {
      const int64_t rowIndex0 = index(x, y, rB0);
      const int64_t rowIndex1 = index(x, y, rB1);

      // block 00
      triplets.emplace_back(rowIndex0, index(x, y, cB0), kroneckerDelta);
      // block 01
      triplets.emplace_back(rowIndex0, index(x, y, cB1), weight);
      // block 10
      triplets.emplace_back(rowIndex1, index(x, y, cB0),
                              weightLaplace   * Laplace<Method>::factor[0]
                            + weightPotential * potential.pixel(x, y)
                            + weightEnergy);

      // Add the non-central elements of the Laplace operator
      for (size_t i = 1; i <= Method::margin; ++i)
      {
        const Real factor = weightLaplace * Laplace<Method>::factor[i];
        triplets.emplace_back(rowIndex1, indexS(x + i, y    , cB0), factor);
        triplets.emplace_back(rowIndex1, indexS(x - i, y    , cB0), factor);
        triplets.emplace_back(rowIndex1, indexS(x    , y + i, cB0), factor);
        triplets.emplace_back(rowIndex1, indexS(x    , y - i, cB0), factor);
      }

      // block 11
      triplets.emplace_back(rowIndex1, index(x, y, cB1), kroneckerDelta);
    }
  }
}

template<class Pde>
PropagatorRungeKutta2<Pde>::PropagatorRungeKutta2(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPotential(wave, params, sample),
  mA(mParams.cols(), mParams.rows()),
  mATemp(mParams.cols(), mParams.rows()),
  mAP0(mParams.cols(), mParams.rows()),
  mAPA(mParams.cols(), mParams.rows()),
  mAPB(mParams.cols(), mParams.rows()),
  mAPC(mParams.cols(), mParams.rows()),
  mB(mParams.cols(), mParams.rows()),
  mBTemp(mParams.cols(), mParams.rows()),
  mBP0(mParams.cols(), mParams.rows()),
  mBPA(mParams.cols(), mParams.rows()),
  mBPB(mParams.cols(), mParams.rows()),
  mBPC(mParams.cols(), mParams.rows()),
  mBandlimitA(CBandlimit::create(mParams.bandlimit(), mParams, mWave, mWave)),
  mBandlimitB(CBandlimit::create(mParams.bandlimit(), mParams, mB, mB)),
  mPde(mParams)
{
  mPde.initB(mA, mB, mPotential);
}

template<class Pde>
void PropagatorRungeKutta2<Pde>::propagate()
{
  // copy the incoming wave function in a buffer with a margin
  mA.assign(mWave);

  const size_t numSlices = mParams.numSlices();
  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    const Real deltaZ = mParams.zEnd(i) - mParams.zStart(i);

    // the actual classical Runge-Kutta algorithm
    mAP0.assign(mB);
    mPde.f(mA, mB, mPotential->potential(2 * i), mBP0);
    mATemp.assign(mA + R(0.5) * deltaZ * mAP0);
    mBTemp.assign(mB + R(0.5) * deltaZ * mBP0);

    mAPA.assign(mBTemp);
    mPde.f(mATemp, mBTemp, mPotential->potential(2 * i + 1), mBPA);
    mATemp.assign(mA + R(0.5) * deltaZ * mAPA);
    mBTemp.assign(mB + R(0.5) * deltaZ * mBPA);

    mAPB.assign(mBTemp);
    mPde.f(mATemp, mBTemp, mPotential->potential(2 * i + 1), mBPB);
    mATemp.assign(mA + deltaZ * mAPB);
    mBTemp.assign(mB + deltaZ * mBPB);

    mAPC.assign(mBTemp);
    mPde.f(mATemp, mBTemp, mPotential->potential(2 * i + 2), mBPC);

    // Combine all results
    mA += deltaZ / R(6.0) * (mAP0 + R(2.0) * (mAPA + mAPB) + mAPC);
    mB += deltaZ / R(6.0) * (mBP0 + R(2.0) * (mBPA + mBPB) + mBPC);

    if (Bandlimit::Disabled != mParams.bandlimit())
    {
      // implement the band-limit
      // copy mA to mWave, as this buffer has no margin and hence can be Fourier
      // transformed
      mWave.assign(mA);
      mBandlimitA->apply();
      mA.assign(mWave);

      mBandlimitB->apply();
    }

    if (verbosity >= 2)
    {
      mPde.diagnostic(i, mZ, mA, mB);
      mZ += deltaZ;
    }

    mWave.assign(mA);

    if (mParams.debug(Debug::propagator))
    {
      static int index = 0;
      if ((index < 200) || (index % 200 == 0))
      {
        mWave.save(mParams.debugDir()  + "/A" + std::to_string(index) + ".tif",
                   mParams, "");
        mB.save(mParams.debugDir()  + "/B" + std::to_string(index) + ".tif",
                mParams, "");
      }
      ++index;
    }

    onWave(i, numSlices, mWave);
  }
}

PropagatorRungeKutta2FT::PropagatorRungeKutta2FT(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPotential(wave, params, sample),
  mA(wave),
  mATemp(mParams.cols(), mParams.rows()),
  mAPA(mParams.cols(), mParams.rows()),
  mAPB(mParams.cols(), mParams.rows()),
  mAPC(mParams.cols(), mParams.rows()),
  mB(mParams.cols(), mParams.rows()),
  mBP0(mParams.cols(), mParams.rows()),
  mBPA(mParams.cols(), mParams.rows()),
  mBPB(mParams.cols(), mParams.rows()),
  mBPC(mParams.cols(), mParams.rows()),
  mTemp(mParams.cols(), mParams.rows()),
  mLaplaceFT(mParams, mTemp, mTemp, Complex(1.0)),
  mWaveNum(mParams.wavenumber()),
  mPotentialFactor(2 * Phys::mE * mParams.lorentzFactor() /
                   powerOf<2>(Phys::hBar)),
  mBandlimitA(CBandlimit::create(mParams.bandlimit(), mParams, mA, mA)),
  mBandlimitB(CBandlimit::create(mParams.bandlimit(), mParams, mB, mB))
{
  mB.setZero();
}

void PropagatorRungeKutta2FT::f(const CBuffer2D& inA, const CBuffer2D& inB,
                                const RBuffer2D& potential, CBuffer2D& outBP)
{
  mTemp.assign(inA);
  mLaplaceFT.apply();
  outBP.assign(- mWaveNum * R(2.0i) * inB
               + mTemp
               + mPotentialFactor * potential * inA);
}

void PropagatorRungeKutta2FT::diagnostic(size_t slice, Real z,
                                         const CBuffer2D& a,
                                         const CBuffer2D& b) const
{
  std::cout << "Slice " << slice << ": " << "z = " << z
         << "; mA.absSqr() = " << a.absSqrReduce() / (mParams.numPixels<Real>())
         << "; mB.absSqr() = " << b.absSqrReduce() / (mParams.numPixels<Real>())
         << '\n';
}

void PropagatorRungeKutta2FT::propagate()
{
  const size_t numSlices = mParams.numSlices();
  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    const Real deltaZ = mParams.zEnd(i) - mParams.zStart(i);

    // the actual classical Runge-Kutta algorithm
    f(mA, mB, mPotential->potential(2 * i), mBP0);
    mATemp.assign(mA + R(0.5) * deltaZ * mB);
    // this used to be mBTemp
    mAPA.assign(mB + R(0.5) * deltaZ * mBP0);

    f(mATemp, mAPA, mPotential->potential(2 * i + 1), mBPA);
    mATemp.assign(mA + R(0.5) * deltaZ * mAPA);
    // this used to be mBTemp
    mAPB.assign(mB + R(0.5) * deltaZ * mBPA);

    f(mATemp, mAPB, mPotential->potential(2 * i + 1), mBPB);
    mATemp.assign(mA + deltaZ * mAPB);
    // this used to be mBTemp
    mAPC.assign(mB + deltaZ * mBPB);

    f(mATemp, mAPC, mPotential->potential(2 * i + 2), mBPC);

    // Combine all results
    mA += deltaZ / 6 * (mB + 2 * (mAPA + mAPB) + mAPC);
    mB += deltaZ / 6 * (mBP0 + 2 * (mBPA + mBPB) + mBPC);

    if (Bandlimit::Disabled != mParams.bandlimit())
    {
      // implement the band-limit
      mBandlimitA->apply();
      mBandlimitB->apply();
    }

    if (verbosity >= 2)
    {
      diagnostic(i, mZ, mA, mB);
      mZ += deltaZ;
    }

    if (mParams.debug(Debug::propagator))
    {
      static int index = 0;
      if ((index < 200) || (index % 200 == 0))
      {
        mA.save(mParams.debugDir() + "/A" + std::to_string(index) + ".tif",
                mParams, "");
        mB.save(mParams.debugDir() + "/B" + std::to_string(index) + ".tif",
                mParams, "");
      }
      ++index;
    }

    onWave(i, numSlices, mA);
  }
}

inline void PdeBase::bufferToEigenVector(const CBuffer2D& buffer,
                                         Vector& vector, int64_t part)
{
  #pragma omp parallel for
  for (int64_t y = 0; y < mRows; ++y)
    for (int64_t x = 0; x < mCols; ++x)
      vector(index(x, y, part)) = buffer.pixel(x, y);
}

inline void PdeBase::eigenVectorToBuffer(const Vector& vector,
                                         CBuffer2D& buffer, int64_t part)
{
  #pragma omp parallel for
  for (int64_t y = 0; y < mRows; ++y)
    for (int64_t x = 0; x < mCols; ++x)
      buffer.pixel(x, y) = Complex(vector(index(x, y, part)));
}

template<class Pde>
PropagatorImplicit1<Pde>::PropagatorImplicit1(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPotential(wave, params, sample),
  mA(mParams.cols(), mParams.rows()),
  mTemp(wave.createCompatibleBuffer()),
  mBandlimitWave(CBandlimit::create(mParams.bandlimit(), mParams, wave, wave)),
  mPde(mParams)
{}

template<class Pde>
void PropagatorTrapezoidal1<Pde>::propagate()
{
  const int64_t numPixels = this->mParams.numPixels();

  // matrix * vectorX = vectorB
  PdeBase::SparseMatrix matrix(numPixels, numPixels);
  PdeBase::Vector vectorB(numPixels);
  PdeBase::Vector vectorX(numPixels);

  PdeBase::Triplets triplets;

  const size_t numSlices = this->mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    const Real deltaZ = this->mParams.zEnd(i) - this->mParams.zStart(i);

    auto& pot0 = this->mPotential->potential(2 * i + 0);
    auto& pot2 = this->mPotential->potential(2 * i + 2);

    // construct the right hand side (vectorB)
    // Copy the wave into mA
    this->mA.assign(this->mWave);

    this->mPde.g(this->mA, pot0, R(0.5) * deltaZ, this->mTemp);
    // copy into vectorB
    this->mPde.bufferToEigenVector(this->mTemp, vectorB);

    // the right hand side is now constructed

    // setting up the matrix
    triplets.clear();
    this->mPde.fillTriplets(triplets, 0, 0, R(-0.5) * deltaZ, pot2);

    matrix.setFromTriplets(triplets.begin(), triplets.end());

    // solve the equation
    vectorX = Eigen::BiCGSTAB<PdeBase::SparseMatrix>(matrix).solve(vectorB);

    // copy back
    this->mPde.eigenVectorToBuffer(vectorX, this->mWave);

    // implement the band-limit
    this->mBandlimitWave->apply();

    if (verbosity >= 2)
    {
      this->mPde.diagnostic(i, this->mZ, this->mA);
      this->mZ += deltaZ;
    }

    this->onWave(i, numSlices, this->mWave);
  }
}

template<class Pde>
void PropagatorLobatto1<Pde>::propagate()
{
  const int64_t cols = this->mParams.cols();
  const int64_t rows = this->mParams.rows();
  const int64_t numPixels = this->mParams.numPixels();

  // matrix * vectorX = vectorB
  PdeBase::SparseMatrix matrix(3 * numPixels, 3 * numPixels);
  PdeBase::Vector vectorX(3 * numPixels);
  PdeBase::Vector vectorB(3 * numPixels);

  PdeBase::Triplets triplets;

  const size_t numSlices = this->mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    const Real deltaZ = this->mParams.zEnd(i) - this->mParams.zStart(i);

    auto& pot0 = this->mPotential->potential(2 * i + 0);
    auto& pot1 = this->mPotential->potential(2 * i + 1);
    auto& pot2 = this->mPotential->potential(2 * i + 2);

    // construct the right hand side
    // Copy the wave into mA
    this->mA.assign(this->mWave);

    // first part of vectorB
    this->mPde.f(this->mA, pot0, this->mTemp);
    this->mPde.bufferToEigenVector(this->mTemp, vectorB, 0);

    // second part of vectorB
    this->mPde.f(this->mA, pot1, this->mTemp);
    this->mPde.bufferToEigenVector(this->mTemp, vectorB, 1);

    // third part of vectorC
    this->mPde.f(this->mA, pot2, this->mTemp);
    this->mPde.bufferToEigenVector(this->mTemp, vectorB, 2);

    // the right hand side is now constructed

    // set up the matrix
    triplets.clear();
#if 0
    // Lobatto IIIc
    this->mPde.fillTriplets(triplets, 0, 0, -1.0 /  6.0 * deltaZ, pot0);
    this->mPde.fillTriplets(triplets, 0, 1,  1.0 /  3.0 * deltaZ, pot0);
    this->mPde.fillTriplets(triplets, 0, 2, -1.0 /  6.0 * deltaZ, pot0);
    this->mPde.fillTriplets(triplets, 1, 0, -1.0 /  6.0 * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 1, 1, -5.0 / 12.0 * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 1, 2,  1.0 / 12.0 * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 2, 0, -1.0 /  6.0 * deltaZ, pot2);
    this->mPde.fillTriplets(triplets, 2, 1, -2.0 /  3.0 * deltaZ, pot2);
    this->mPde.fillTriplets(triplets, 2, 2, -1.0 /  6.0 * deltaZ, pot2);
#else
    // Lobatto IIIa
    this->mPde.fillTriplets(triplets, 0, 0, R( 0.0)           * deltaZ, pot0);
    this->mPde.fillTriplets(triplets, 1, 0, R(-5.0) / R(24.0) * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 1, 1, R(-1.0) / R( 3.0) * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 1, 2, R( 1.0) / R(24.0) * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 2, 0, R(-1.0) / R( 6.0) * deltaZ, pot2);
    this->mPde.fillTriplets(triplets, 2, 1, R(-2.0) / R( 3.0) * deltaZ, pot2);
    this->mPde.fillTriplets(triplets, 2, 2, R(-1.0) / R( 6.0) * deltaZ, pot2);
#endif

    matrix.setFromTriplets(triplets.begin(), triplets.end());

    // solve the equation
    vectorX = Eigen::BiCGSTAB<PdeBase::SparseMatrix>(matrix).solve(vectorB);

    // calculate the final result in store it in mWave
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
    {
      for (int64_t x = 0; x < cols; ++x)
      {
        auto sum =          vectorX(this->mPde.index(x, y, 0)) +
                   R(4.0) * vectorX(this->mPde.index(x, y, 1)) +
                            vectorX(this->mPde.index(x, y, 2));

        this->mWave.pixel(x, y) += deltaZ / R(6.0) * Complex(sum);
      }
    }

    // implement the band-limit
    this->mBandlimitWave->apply();

    if (verbosity >= 2)
    {
      this->mPde.diagnostic(i, this->mZ, this->mA);
      this->mZ += deltaZ;
    }

    this->onWave(i, numSlices, this->mWave);
  }
}

template<class Pde>
PropagatorImplicit2<Pde>::PropagatorImplicit2(CBuffer2D& wave,
                                             const MultisliceParameters& params,
                                             const ConstSamplePtr& sample) :
  PropagatorPotential(wave, params, sample),
  mA(mParams.cols(), mParams.rows()),
  mB(mWave.createCompatibleBuffer()),
  mTemp(wave.createCompatibleBuffer()),
  mBandlimitWave(CBandlimit::create(mParams.bandlimit(), mParams, mWave, mWave)),
  mBandlimitB(CBandlimit::create(mParams.bandlimit(), mParams, mB, mB)),
  mPde(mParams)
{
  mPde.initB(mA, mB, mPotential);
}

template<class Pde>
void PropagatorTrapezoidal2<Pde>::propagate()
{
  const int64_t numPixels = this->mParams.numPixels();

  // matrix * vectorX = vectorB
  PdeBase::SparseMatrix matrix(2 * numPixels, 2 * numPixels);
  PdeBase::Vector vectorX(2 * numPixels);
  PdeBase::Vector vectorB(2 * numPixels);

  PdeBase::Triplets triplets;

  const size_t numSlices = this->mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    const Real deltaZ = this->mParams.zEnd(i) - this->mParams.zStart(i);

    auto& pot0 = this->mPotential->potential(2 * i + 0);
    auto& pot2 = this->mPotential->potential(2 * i + 2);

    // construct the right hand side (vectorB)
    // Copy the wave into mA
    this->mA.assign(this->mWave);

    // fill the first part of vectorB
    this->mTemp.assign(this->mA + R(0.5) * deltaZ * this->mB);
    this->mPde.bufferToEigenVector(this->mTemp, vectorB, 0);

    // fill the second part of vectorB
    this->mPde.g(this->mA, this->mB, pot0, R(0.5) * deltaZ, this->mTemp);
    this->mPde.bufferToEigenVector(this->mTemp, vectorB, 1);

    // the right hand side is now constructed

    // setting up the mAtrix
    triplets.clear();
    this->mPde.fillTriplets(triplets, 0, 0, R(-0.5) * deltaZ, pot2);

    matrix.setFromTriplets(triplets.begin(), triplets.end());

    // solve the equation
    vectorX = Eigen::BiCGSTAB<PdeBase::SparseMatrix>(matrix).solve(vectorB);

    // copy the first part to mWave
    this->mPde.eigenVectorToBuffer(vectorX, this->mWave, 0);

    // copy the second part to mB
    this->mPde.eigenVectorToBuffer(vectorX, this->mB, 1);

    // implement the band-limit
    this->mBandlimitWave->apply();
    this->mBandlimitB->apply();

    if (verbosity >= 2)
    {
      this->mPde.diagnostic(i, this->mZ, this->mA, this->mB);
      this->mZ += deltaZ;
    }

    this->onWave(i, numSlices, this->mWave);
  }
}

template<class Pde>
void PropagatorLobatto2<Pde>::propagate()
{
  const int64_t cols = this->mParams.cols();
  const int64_t rows = this->mParams.rows();
  const int64_t numPixels = this->mParams.numPixels();

  // matrix * vectorX = vectorB
  PdeBase::SparseMatrix matrix(6 * numPixels, 6 * numPixels);
  PdeBase::Vector vectorX(6 * numPixels);
  PdeBase::Vector vectorB(6 * numPixels);

  PdeBase::Triplets triplets;

  const size_t numSlices = this->mParams.numSlices();

  // for all slices
  for (size_t i = 0; i < numSlices; ++i)
  {
    const Real deltaZ = this->mParams.zEnd(i) - this->mParams.zStart(i);

    auto& pot0 = this->mPotential->potential(2 * i + 0);
    auto& pot1 = this->mPotential->potential(2 * i + 1);
    auto& pot2 = this->mPotential->potential(2 * i + 2);

    // construct the right hand side (vectorB)
    // Copy the wave into mA
    this->mA.assign(this->mWave);

    // Copy mB into the 1st, 3rd, 5th part of vectorB
    this->mPde.bufferToEigenVector(this->mB, vectorB, 0);
    this->mPde.bufferToEigenVector(this->mB, vectorB, 2);
    this->mPde.bufferToEigenVector(this->mB, vectorB, 4);

    // fill the second part of vectorB
    this->mPde.f(this->mA, this->mB, pot0, this->mTemp);
    this->mPde.bufferToEigenVector(this->mTemp, vectorB, 1);

    // fill the fourth part of vectorB
    this->mPde.f(this->mA, this->mB, pot1, this->mTemp);
    this->mPde.bufferToEigenVector(this->mTemp, vectorB, 3);

    // fill the sixth part of vectorB
    this->mPde.f(this->mA, this->mB, pot0, this->mTemp);
    this->mPde.bufferToEigenVector(this->mTemp, vectorB, 5);

    // the right hand side is now constructed

    // set up the matrix
    triplets.clear();
#if 0
    // Lobatto IIIc
    this->mPde.fillTriplets(triplets, 0, 0, -1.0 /  6.0 * deltaZ, pot0);
    this->mPde.fillTriplets(triplets, 0, 1,  1.0 /  3.0 * deltaZ, pot0);
    this->mPde.fillTriplets(triplets, 0, 2, -1.0 /  6.0 * deltaZ, pot0);
    this->mPde.fillTriplets(triplets, 1, 0, -1.0 /  6.0 * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 1, 1, -5.0 / 12.0 * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 1, 2,  1.0 / 12.0 * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 2, 0, -1.0 /  6.0 * deltaZ, pot2);
    this->mPde.fillTriplets(triplets, 2, 1, -2.0 /  3.0 * deltaZ, pot2);
    this->mPde.fillTriplets(triplets, 2, 2, -1.0 /  6.0 * deltaZ, pot2);
#else
    // Lobatto IIIa
    this->mPde.fillTriplets(triplets, 0, 0, R( 0.0)           * deltaZ, pot0);
    this->mPde.fillTriplets(triplets, 1, 0, R(-5.0) / R(24.0) * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 1, 1, R(-1.0) / R( 3.0) * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 1, 2, R( 1.0) / R(24.0) * deltaZ, pot1);
    this->mPde.fillTriplets(triplets, 2, 0, R(-1.0) / R( 6.0) * deltaZ, pot2);
    this->mPde.fillTriplets(triplets, 2, 1, R(-2.0) / R( 3.0) * deltaZ, pot2);
    this->mPde.fillTriplets(triplets, 2, 2, R(-1.0) / R( 6.0) * deltaZ, pot2);
#endif

    matrix.setFromTriplets(triplets.begin(), triplets.end());

    // solve the equation
    vectorX = Eigen::BiCGSTAB<PdeBase::SparseMatrix>(matrix).solve(vectorB);

    // calculate the final result and store it in mWave
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
    {
      for (int64_t x = 0; x < cols; ++x)
      {
        auto sum =          vectorX(this->mPde.index(x, y, 0)) +
                   R(4.0) * vectorX(this->mPde.index(x, y, 2)) +
                            vectorX(this->mPde.index(x, y, 4));

        this->mWave.pixel(x, y) += deltaZ / 6 * Complex(sum);
      }
    }

    // the same for mB
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
    {
      for (int64_t x = 0; x < cols; ++x)
      {
        auto sum =          vectorX(this->mPde.index(x, y, 1)) +
                   R(4.0) * vectorX(this->mPde.index(x, y, 3)) +
                            vectorX(this->mPde.index(x, y, 5));

        this->mB.pixel(x, y) += deltaZ / 6 * Complex(sum);
      }
    }

    // implement the band-limit
    this->mBandlimitWave->apply();
    this->mBandlimitB->apply();

    if (verbosity >= 2)
    {
      this->mPde.diagnostic(i, this->mZ, this->mA, this->mB);
      this->mZ += deltaZ;
    }

    if (this->mParams.debug(Debug::propagator))
    {
      static int index = 0;
      if ((index < 200) || (index % 200 == 0))
      {
        std::string debugDir = this->mParams.debugDir();
        this->mWave.save(debugDir + "/A" + std::to_string(index) + ".tiff",
                         this->mParams, "");
        this->mB.save(debugDir + "/B" + std::to_string(index) + ".tiff",
                      this->mParams, "");
      }
      ++index;
    }

    this->onWave(i, numSlices, this->mWave);
  }
}

class RegisterPropagators
{
public:
  RegisterPropagators()
  {
    Propagator::registerClass<PropagatorSingle<ThreePoint>>("Single-3Pt");
    Propagator::registerClass<PropagatorSingle<FivePoint>> ("Single-5Pt");
    Propagator::registerClass<PropagatorSingle<SevenPoint>>("Single-7Pt");
    Propagator::registerClass<PropagatorSingle<NinePoint>> ("Single-9Pt");
    Propagator::registerClass<PropagatorSingleFT>          ("Single-FT");
    Propagator::registerClass<PropagatorClassical<ThreePoint>>("Classical-3Pt");
    Propagator::registerClass<PropagatorClassical<FivePoint>> ("Classical-5Pt");
    Propagator::registerClass<PropagatorClassical<SevenPoint>>("Classical-7Pt");
    Propagator::registerClass<PropagatorClassical<NinePoint>> ("Classical-9Pt");
    Propagator::registerClass<PropagatorClassicalFT>          ("Classical-FT");
    Propagator::registerClass<PropagatorPartial<ThreePoint>>("Partial-3Pt");
    Propagator::registerClass<PropagatorPartial<FivePoint>> ("Partial-5Pt");
    Propagator::registerClass<PropagatorPartial<SevenPoint>>("Partial-7Pt");
    Propagator::registerClass<PropagatorPartial<NinePoint>> ("Partial-9Pt");
    Propagator::registerClass<PropagatorFull<ThreePoint>>("Full-3Pt");
    Propagator::registerClass<PropagatorFull<FivePoint>> ("Full-5Pt");
    Propagator::registerClass<PropagatorFull<SevenPoint>>("Full-7Pt");
    Propagator::registerClass<PropagatorFull<NinePoint>> ("Full-9Pt");
    Propagator::registerClass<PropagatorFullFT>          ("Full-FT");
    Propagator::registerClass<PropagatorFullFT2>         ("Full-FT2");
    Propagator::registerClass<PropagatorRungeKutta1<Pde1<ThreePoint>>>("RungeKutta1-3Pt");
    Propagator::registerClass<PropagatorRungeKutta1<Pde1<FivePoint>>> ("RungeKutta1-5Pt");
    Propagator::registerClass<PropagatorRungeKutta1<Pde1<SevenPoint>>>("RungeKutta1-7Pt");
    Propagator::registerClass<PropagatorRungeKutta1<Pde1<NinePoint>>> ("RungeKutta1-9Pt");
    Propagator::registerClass<PropagatorRungeKutta1FT>                ("RungeKutta1-FT");
    Propagator::registerClass<PropagatorRungeKutta2<Pde2<ThreePoint>>>("RungeKutta2-3Pt");
    Propagator::registerClass<PropagatorRungeKutta2<Pde2<FivePoint>>> ("RungeKutta2-5Pt");
    Propagator::registerClass<PropagatorRungeKutta2<Pde2<SevenPoint>>>("RungeKutta2-7Pt");
    Propagator::registerClass<PropagatorRungeKutta2<Pde2<NinePoint>>> ("RungeKutta2-9Pt");
    Propagator::registerClass<PropagatorRungeKutta2FT>                ("RungeKutta2-FT");
    Propagator::registerClass<PropagatorRungeKutta2<Pde3<NinePoint>>> ("RungeKutta3-9Pt");
    Propagator::registerClass<PropagatorTrapezoidal1<Pde1<ThreePoint>>>("Trapezoidal1-3Pt");
    Propagator::registerClass<PropagatorTrapezoidal1<Pde1<FivePoint>>> ("Trapezoidal1-5Pt");
    Propagator::registerClass<PropagatorTrapezoidal1<Pde1<SevenPoint>>>("Trapezoidal1-7Pt");
    Propagator::registerClass<PropagatorTrapezoidal1<Pde1<NinePoint>>> ("Trapezoidal1-9Pt");
    Propagator::registerClass<PropagatorTrapezoidal2<Pde2<ThreePoint>>>("Trapezoidal2-3Pt");
    Propagator::registerClass<PropagatorTrapezoidal2<Pde2<FivePoint>>> ("Trapezoidal2-5Pt");
    Propagator::registerClass<PropagatorTrapezoidal2<Pde2<SevenPoint>>>("Trapezoidal2-7Pt");
    Propagator::registerClass<PropagatorTrapezoidal2<Pde2<NinePoint>>> ("Trapezoidal2-9Pt");
    Propagator::registerClass<PropagatorLobatto1<Pde1<ThreePoint>>>("Lobatto1-3Pt");
    Propagator::registerClass<PropagatorLobatto1<Pde1<FivePoint>>> ("Lobatto1-5Pt");
    Propagator::registerClass<PropagatorLobatto1<Pde1<ThreePoint>>>("Lobatto1-7Pt");
    Propagator::registerClass<PropagatorLobatto1<Pde1<NinePoint>>> ("Lobatto1-9Pt");
    Propagator::registerClass<PropagatorLobatto2<Pde2<ThreePoint>>>("Lobatto2-3Pt");
    Propagator::registerClass<PropagatorLobatto2<Pde2<FivePoint>>> ("Lobatto2-5Pt");
    Propagator::registerClass<PropagatorLobatto2<Pde2<SevenPoint>>>("Lobatto2-7Pt");
    Propagator::registerClass<PropagatorLobatto2<Pde2<NinePoint>>> ("Lobatto2-9Pt");
  }
} registerPropagators;

} // namespace Aurora
