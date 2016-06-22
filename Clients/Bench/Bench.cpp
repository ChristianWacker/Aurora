//--- Aurora/Clients/Bench/Bench.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Crystal.hpp"
#include "AuroraLib/Ctf.hpp"
#include "AuroraLib/Fft.hpp"

#if (AURORA_USE_FFTW == 1)
# include "AuroraLib/FftFftw.hpp"
#endif

#include "AuroraLib/FftKiss.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Lapack.hpp"
#include "AuroraLib/MatrixSupport.hpp"
#include "AuroraLib/MultisliceParameters.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/Statistics.hpp"
#include "AuroraLib/TemSimulation.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace Aurora
{

enum class Mode
{
  fft, tem, lu, eigenvalue, matMul
};

CL::Option<Mode> mode("mode", "Benchmark mode", Mode::tem,
  enumVal(Mode::fft), enumVal(Mode::tem), enumVal(Mode::lu),
  enumVal(Mode::eigenvalue), enumVal(Mode::matMul));
CL::Option<int> dimension("dimension", "", 1000);

Real gauss(Real x, Real sigma, Real mu)
{
  return std::exp(-R(0.5) * powerOf<2>(x - mu) / sigma) /
    (sigma * std::sqrt(Math::twoPi));
}

/// Base class for all benchmarks
class Benchmark
{
public:
  virtual void bench()
  {
    std::cout << "Benchmarking...\n";

    Statistics<std::chrono::high_resolution_clock::duration> statistics;

    for (size_t i = 0; i < 25; ++i)
    {
      setup();

      auto start = std::chrono::high_resolution_clock::now();
      run();
      auto stop = std::chrono::high_resolution_clock::now();

      teardown();

      if (i >= 10)
        statistics.add(stop - start);
    }

    std::cout << "Time needed - Mean: "
              << duration(statistics.mean())
              << " Median: " << duration(statistics.median()) << std::endl;
  }

  virtual ~Benchmark() {}

private:
  virtual void run() {}
  virtual void setup() {}
  virtual void teardown() {}
};

class BenchmarkFft : public Benchmark
{
  template<class T>
  void benchFft(bool skipNonPowerOf2 = false)
  {
    const size_t repetitions = 25;

    // test different resolutions
    for (int64_t resolution : {256, 384, 512, 741, 768, 1024, 2048})
    {
      // include chrono for time measuring
      using namespace std::chrono;

      if (skipNonPowerOf2 && !isPowerOfTwo(resolution))
        continue;

      std::cout << "- " << resolution << "x" << resolution << " - ";

      // acquire memory
      CBuffer2D bufferIn(resolution, resolution);
      CBuffer2D bufferOut(resolution, resolution);

      // set up FFT
      auto startInit = high_resolution_clock::now();
      auto fft = std::make_unique<T>(bufferIn, bufferOut, Direction::forward);
      auto stopInit = high_resolution_clock::now();

      std::cout << "Time needed for init: "
                << duration(stopInit - startInit) << " - ";

      // generate test data
      const Real sigma = 1.0;
      for (int64_t y = 0; y < resolution; ++y)
      {
        for (int64_t x = 0; x < resolution; ++x)
        {
          // normalized coordinates [-1.0, 1.0]
          const Real xNorm = 2 * ((Real(x) / Real(resolution - 1)) - R(0.5));
          const Real yNorm = 2 * ((Real(y) / Real(resolution - 1)) - R(0.5));

          Real distanceToCenter = std::sqrt(xNorm * xNorm + yNorm * yNorm);

          bufferIn.pixel(x, y) = gauss(distanceToCenter, sigma, 0);
        }
      }

      Statistics<high_resolution_clock::duration> runStatistics;

      for (size_t i = 0; i < repetitions + 10; ++i)
      {
        auto startTransform = high_resolution_clock::now();
        fft->transform();
        auto stopTransform = high_resolution_clock::now();

        if (i >= 10)
          runStatistics.add(stopTransform - startTransform);
      }

      std::cout << "Time needed for FFT -"
                   " Mean: " << duration(runStatistics.mean())
                << " Median: " << duration(runStatistics.median())
                << std::endl;
    }
  }

  void bench() override
  {
    std::cout << "Benchmarking two dimensional FFT...\n";

    std::cout << "Bench Kiss-FFT:\n";
    benchFft<FftKiss>();

  #if (1 == AURORA_USE_FFTW)
    std::cout << "Bench FFTW-Estimate:\n";
    benchFft<FftFftw<FFTW_ESTIMATE>>();
    std::cout << "Bench FFTW-Measure:\n";
    benchFft<FftFftw<FFTW_MEASURE>>();
    std::cout << "Bench FFTW-Patient:\n";
    benchFft<FftFftw<FFTW_PATIENT>>();
  #endif
  }
};

class BenchmarkTemSimulation : public Benchmark
{
  void bench() override
  {
    std::cout << "Benchmarking TEM simulation...\n";

    // suppress output
    verbosity.setValue(0);

    // create a silicon sample
    auto siliconCrystal = Crystal::create("Silicon");

    // rotate the sample in [110] direction
    Matrix3x3R rotation;
    rotation << R(0.0), R(0.0), R(1.0),
                R(1.0) / std::sqrt(R(2.0)), R(-1.0) / std::sqrt(R(2.0)), R(0.0),
                R(1.0) / std::sqrt(R(2.0)), R( 1.0) / std::sqrt(R(2.0)), R(0.0);
    auto siliconSample =
      siliconCrystal->generateSample(Vector3R(-R(10.0), -R(10.0), -R(10.0)),
                                     Vector3R( R(10.0),  R(10.0),  R(10.0)),
                                     rotation,
                                     Vector3R::Zero());

    siliconSample->preprocess();

    Statistics<std::chrono::high_resolution_clock::duration> statistics;
    for (size_t i = 0; i < 15; ++i)
    {
      auto start = std::chrono::high_resolution_clock::now();

      MultisliceParameters params;
      params.setResolution(512, 512);
      params.setEnergy(R(1e5));
      params.setBox(Vector3R(-R(5.43), -R(3.84), -R(4.0)),
                    Vector3R( R(5.43),  R(3.84),  R(4.0)),
                    R(1.0));
      params.setCutoff(R(5.0));
      params.setFormFactorName("Peng");
      params.setPhaseShiftName("Opt");
      params.setBandlimit(Bandlimit::Smooth);

      TemSimulation tem(params, 0, siliconSample, "Classical-FT");

      using namespace std::placeholders;
      tem.onWave = std::bind(&BenchmarkTemSimulation::onWave, this, _1, _2, _3);

      tem.run();

      auto stop = std::chrono::high_resolution_clock::now();

      if (i >= 5)
        statistics.add(stop - start);
    }

    std::cout << "Time needed for TEM simulation - Mean: "
              << duration(statistics.mean())
              << " Median: " << duration(statistics.median()) << std::endl;
  }

  void onWave(size_t slice, size_t numSlices, const CBuffer2D& wave)
  {
    if (numSlices - 1 == slice)
    {
      const int64_t cols = wave.cols();
      const int64_t rows = wave.rows();

      auto waveFft = wave.fft(Direction::forward);

      CtfCoherent ctf(1e5);
      auto ctfBuffer = ctf.buffer(Buffer2DInfoBase(cols, rows, R(10.86), R(7.68)),
                                  R(1.0) / (cols * rows));
      auto waveAfterOpticsFft = waveFft.clone();
      waveAfterOpticsFft *= ctfBuffer;

      auto waveAfterOptics(waveAfterOpticsFft.fft(Direction::backward));

      // calculate the total probability
      const Real norm = waveAfterOptics.absSqrReduce();

      std::cout << "Probability = " << norm / (cols * rows) << '\n';
    }
  }
};

class BenchmarkEVDecomposition : public Benchmark
{
  void setup() override
  {
    mMat = MatrixNC::Random(dimension, dimension);
    mMat += mMat.adjoint().eval();
  }

  void run() override
  {
    SelfAdjointEigenSolver ev(mMat);
    std::cout << ev.eigenvalues().norm() + ev.eigenvectors().norm() << '\n';
  }

  MatrixNC mMat;
};

class BenchmarkLUDecomposition : public Benchmark
{
  void setup() override
  {
    mMat = MatrixNC::Random(dimension, dimension);
  }

  void run() override
  {
    PartialPivLU lu(mMat);
    std::cout << lu.matrixLU().norm() << '\n';
  }

  MatrixNC mMat;
};

class BenchmarkMatrixMultiplication : public Benchmark
{
  void setup() override
  {
    mA = MatrixNC::Random(dimension, dimension);
    mB = MatrixNC::Random(dimension, dimension);
  }

  void run() override
  {
    MatrixNC c = mA * mB;
    std::cout << c.norm() << '\n';
  }

  MatrixNC mA, mB;
};

class Application
{
public:
  void benchFft();
  void benchTemSimulation();
  void run(int argc, char** argv);

private:
  void onWave(size_t slice, size_t numSlices, const CBuffer2D& waveFunction);
};

void Application::run(int argc, char** argv)
{
  std::cout << "Benchmark\n";

  if (verbosity >= 1)
    std::cout << buildString() << '\n';

  if (!CL::parse(argc, argv))
    return;

  std::unique_ptr<Benchmark> benchmark;
  switch (mode)
  {
  case Mode::fft:
    benchmark = std::make_unique<BenchmarkFft>();
    break;
  case Mode::tem:
    benchmark = std::make_unique<BenchmarkTemSimulation>();
    break;
  case Mode::eigenvalue:
    benchmark = std::make_unique<BenchmarkEVDecomposition>();
    break;
  case Mode::lu:
    benchmark = std::make_unique<BenchmarkLUDecomposition>();
    break;
  case Mode::matMul:
    benchmark = std::make_unique<BenchmarkMatrixMultiplication>();
    break;
  }
  benchmark->bench();
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Application application;
    application.run(argc, argv);
  }
  catch (Aurora::Exception& e)
  {
    std::cerr << "An exception occured:\n" << e.string() << "\n";
  }
}
