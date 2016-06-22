//--- Aurora/Clients/Matrix/PassPower.cpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "PassPower.hpp"

#include <chrono>
#include <iostream>

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"

namespace Aurora
{

PassPower::PassPower(const MultisliceParameters& params,
    const std::string& matrixDir, size_t power,
    const std::string& matrixPrefix, size_t startIndex, size_t increment) :
  Pass(params, matrixDir), mPower(power), mMatrixPrefix(matrixPrefix),
  mStartIndex(startIndex), mIncrement(increment)
{
  assert(mPower >= 1 && "Power is invalid.");
  assert(mStartIndex <= mPower);
}

void PassPower::run()
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  MatrixNC transferMatrix;
  loadMatrix(matrixName(mMatrixPrefix + "-power", mIncrement), transferMatrix);

  MatrixNC powerMatrix;
  loadMatrix(matrixName(mMatrixPrefix + "-power", mStartIndex), powerMatrix);
  processMatrix(mStartIndex, powerMatrix);

  for (size_t i = mStartIndex + mIncrement; i <= mPower; i += mIncrement)
  {
    powerMatrix *= transferMatrix;
    processMatrix(i, powerMatrix);

    if (verbosity >= 1)
    {
      std::cout  << "Elapsed time (Power " << i << "): "
                 << duration(high_resolution_clock::now() - startTime)
                 << std::endl;
    }
  }
}

void PassPowerF::run()
{
  MatrixNC mat;
  loadMatrix(matrixName("F-power", mPower), mat);

  MatrixNC u;
  loadEigenvectors(matrixName("A-", 0), u);

  // plane wave
  VectorNC f0 = u.adjoint() * VectorNC::Ones(mParams.numPixels());

  loadEigenvectors(matrixName("A-", mParams.numSlices()), u);
  VectorNC exitWave = u * (mat * f0).eval();

  CBuffer2D waveBuffer(mParams.cols(), mParams.rows());
  vectorToBuffer(exitWave, waveBuffer);

  saveWave(waveBuffer, mParams, "Power", mPower);
}

PassPowerT::PassPowerT(const MultisliceParameters& params,
    const std::string& matrixDir, size_t power) :
  PassPower(params, matrixDir, power, "T")
{}

void PassPowerT::processMatrix(size_t power, const MatrixNC& mat)
{
  using namespace std::chrono;
  const int64_t numPix = mParams.numPixels();

  // T11
  saveMatrix(mat.block(0, 0, numPix, numPix),
             matrixDir() + "/T-power" + toString(power, 4) + "-11.mat");

  // T12
  saveMatrix(mat.block(0, numPix, numPix, numPix),
             matrixDir() + "/T-power" + toString(power, 4) + "-12.mat");

  auto startTime = high_resolution_clock::now();

  // todo: optimize LU decomposition and inverse
  // T_22^-1 * T_21
  MatrixNC productMatrix =
    mat.block(numPix, numPix, numPix, numPix).partialPivLu().inverse() *
    mat.block(numPix, 0, numPix, numPix);

  if (verbosity >= 1)
  {
    std::cout << "Time for matrix inversion: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }

  saveMatrix(productMatrix,
             matrixDir() + "/T-power" + toString(power, 4) + "-product.mat");
}

void PassPowerLarge::run()
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  for (size_t i = mStartIndex + mIncrement; i <= mPower; i += mIncrement)
  {
    std::string transferName =
      matrixDir() + "/" + mMatrixPrefix + "-power" + toString(mIncrement, 4);
    std::string powerSrcName =
      matrixDir() + "/" + mMatrixPrefix + "-power" + toString(i - mIncrement, 4);
    std::string powerDstName =
      matrixDir() + "/" + mMatrixPrefix + "-power" + toString(i, 4);

    multiplyLarge(transferName, powerSrcName, powerDstName);

    if (verbosity >= 1)
    {
      std::cout  << "Elapsed time (Repetition " << i << "): "
                 << duration(high_resolution_clock::now() - startTime)
                 << std::endl;
    }
  }
}

} // namespace Aurora
