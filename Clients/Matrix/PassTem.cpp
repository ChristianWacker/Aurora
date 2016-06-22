//--- Aurora/Clients/Matrix/PassTem.cpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "PassTem.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Lapack.hpp"

#include <chrono>

namespace Aurora
{

void PassTem::run()
{
  MatrixNC productMatrix;
  loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) +
             "-product.mat", productMatrix);

  MatrixNC t11;
  loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) + "-11.mat", t11);

  MatrixNC t12;
  loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) + "-12.mat", t12);

  MatrixNC u0;
  loadEigenvectors(matrixName("A-", 0), u0);

  // plane incoming wave
  VectorNC f0 = u0.adjoint() * VectorNC::Ones(mParams.numPixels());
  VectorNC b0 = -productMatrix * f0;

  CBuffer2D temporaryWaveBuffer(mParams.cols(), mParams.rows());
  vectorToBuffer(u0 * b0, temporaryWaveBuffer);
  saveWave(temporaryWaveBuffer, mParams, "backward", mPower);

  MatrixNC uN;
  loadEigenvectors(matrixName("A-", mParams.numSlices()), uN);
  vectorToBuffer(uN * (t11 * f0 + t12 * b0), temporaryWaveBuffer);
  saveWave(temporaryWaveBuffer, mParams, "forward", mPower);
}

void PassTemLarge::run()
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  MatrixNC u0;
  loadEigenvectors(matrixName("A-", 0), u0);

  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();

  // plane incoming wave
  const VectorNC f0 = u0.adjoint() * VectorNC::Ones(cols * rows);

  VectorNC b0;
  {
    PartialPivLU luDecomposition;

    {
      MatrixNC t22;
      loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) +
                 "-22.mat", t22);

      luDecomposition.compute(t22);
    }

    if (verbosity >= 1)
    {
      std::cout << "Time for LU decomposition: "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }

    MatrixNC t21;
    loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) +
               "-21.mat", t21);

    b0 = luDecomposition.solve(-t21 * f0);

    if (verbosity >= 1)
    {
      std::cout << "Time for solve(): "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }

  MatrixNC t11;
  loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) +
             "-11.mat", t11);

  MatrixNC t12;
  loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) +
             "-12.mat", t12);

  CBuffer2D tempWaveBuffer(mParams.cols(), mParams.rows());
  vectorToBuffer(u0 * b0, tempWaveBuffer);
  saveWave(tempWaveBuffer, mParams, "backward", mPower);

  u0.resize(0, 0);
  MatrixNC uN;
  loadEigenvectors(matrixName("A-", mParams.numSlices()), uN);
  vectorToBuffer(uN * ((t11 * f0).eval() + (t12 * b0).eval()), tempWaveBuffer);
  saveWave(tempWaveBuffer, mParams, "forward", mPower);
}

PassTemIterative::PassTemIterative(const MultisliceParameters& params,
  const std::string& matrixDir, size_t power, const std::string& outTemplate,
  size_t minIterations, size_t maxIterations) :
  PassTem(params, matrixDir, power, outTemplate),
  mMinIterations(minIterations), mMaxIterations(maxIterations)
{ }

void PassTemIterative::run()
{
  // load the necessary matrices
  MatrixNC u0;
  loadEigenvectors(matrixName("A-", 0), u0);

  // plane incoming wave
  VectorNC psiIn = VectorNC::Ones(mParams.numPixels());
  VectorNC propagatedF0;

  {
    MatrixNC t11;
    loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) + "-11.mat", t11);
    propagatedF0 = t11 * (u0.adjoint() * psiIn).eval();
  }

  MatrixNC t12;
  loadMatrix(matrixDir() + "/T-power" + toString(mPower, 4) + "-12.mat", t12);

  MatrixNC s21;
  loadMatrix(matrixDir() + "/S-power" + toString(mPower, 4) + "-21.mat", s21);

  // initial values: no backscattered wave
  VectorNC b0Cur = VectorNC::Zero(mParams.numPixels());

  // unmodified incomming wave
  VectorNC fNCur = propagatedF0;

  // the vectors of the previous iteration
  VectorNC b0Prev;
  VectorNC fNPrev;

  size_t i = 0;
  while (true)
  {
    ++i;

    b0Prev = b0Cur;
    fNPrev = fNCur;

    b0Cur = s21 * fNPrev;
    fNCur = t12 * b0Prev + propagatedF0;

    const Real distanceB = (b0Cur - b0Prev).squaredNorm();
    const Real distanceF = (fNCur - fNPrev).squaredNorm();

    if (verbosity >= 1)
    {
      std::cout << "Iteration " << i << ":\n"
                   "- DistanceB = " << distanceB << "\n"
                   "- DistanceF = " << distanceF << '\n';
    }

    if ((i > mMinIterations) && (distanceF < R(1e-9)) && (distanceB < R(1e-9)))
      break;

    if ((i > mMaxIterations) || std::isnan(distanceF) || std::isnan(distanceB))
    {
      std::cout << "Divergent!\n";
      break;
    }
  }

  t12.resize(0, 0);
  s21.resize(0, 0);

  CBuffer2D temporaryWaveBuffer(mParams.cols(), mParams.rows());

  vectorToBuffer(u0 * b0Cur, temporaryWaveBuffer);
  saveWave(temporaryWaveBuffer, mParams, "backward", mPower);

  MatrixNC uN;
  loadEigenvectors(matrixName("A-", mParams.numSlices()), uN);
  vectorToBuffer(uN * fNCur, temporaryWaveBuffer);
  saveWave(temporaryWaveBuffer, mParams, "forward", mPower);
}

} // namespace Aurora
