//--- Aurora/Clients/Matrix/Pass1.cpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "Pass1.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Lapack.hpp"
#include "AuroraLib/PhaseShift.hpp"
#include "AuroraLib/MultisliceParameters.hpp"

#include <chrono>
#include <iostream>

namespace Aurora
{

using namespace std::chrono;

MatrixA2::MatrixA2(const MultisliceParameters &params,
                   const ConstSamplePtr& sample,
                   const std::string& filenameLaplaceMatrix) :
  mParams(params), mFilenameLaplaceMatrix(filenameLaplaceMatrix)
{
  // Use the sample information to calculate the phase shifts
  // Explicitly deactivate the bandlimit, we will do it by ourselves
  mParams.setBandlimit(Bandlimit::Disabled);

  // extend the bounding box slightly so we can also calculate the matrix for
  // next slice after simulation area
  mParams.setBox(mParams.minBox(),
                 mParams.maxBox() + Vector3R(0, 0, mParams.deltaZ()),
                 mParams.deltaZ());
  mPhaseShifts = PhaseShift::create(mParams.phaseShiftName(), mParams, sample,
                                    PhaseShiftMode::phaseShift);
}

void MatrixA2::calc(size_t slice, VectorNR& eigenvalues,
                    MatrixNC& eigenvectors) const
{
  auto startTime = high_resolution_clock::now();
  const int64_t numPixels = mParams.numPixels();

  if (verbosity >= 1)
    std::cout << "   * A (" << slice << ") started..." << std::endl;

  // This is the squared matrix A^2
  MatrixNC a2(numPixels, numPixels);

  {
    calcDiagonalPart(a2, slice);

    // Load the matrix for the Laplace operator
    MatrixNC laplaceMatrix;
    loadMatrix(mFilenameLaplaceMatrix, laplaceMatrix);

    a2 -= laplaceMatrix;
  }

  if (verbosity >= 1)
  {
    std::cout << "   * A (" << slice << ") Laplace matrix applied."
              << std::endl;
  }

  if (mParams.debug(Debug::pass1))
    saveMatrix(a2, mParams.debugDir() + "/A2-" + toString(slice, 4) + ".mat");

  // A^2 is now constructed

  // Calculate the eigen values and eigen vectors
  SelfAdjointEigenSolver eigenSolver(a2);

  if (verbosity >= 1)
  {
    std::cout << "   * A (" << slice << ") Decomposition ";

    if (eigenSolver.info() == Eigen::ComputationInfo::NoConvergence)
      std::cout << "divergent." << std::endl;
    else
      std::cout << "calculated." << std::endl;
  }

  // free the memory
  a2.resize(0, 0);
  eigenvalues = eigenSolver.eigenvalues();
  eigenvectors = eigenSolver.eigenvectors();

  if (verbosity >= 1)
  {
    std::cout << "   * (" << slice << ") Runtime: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void MatrixA2::calcDiagonalPart(Eigen::Ref<MatrixNC> diagonalMatrix,
  size_t sliceIndex) const
{
  auto startTime = high_resolution_clock::now();

  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const int64_t numPixels = cols * rows;
  const Real wavenum = mParams.wavenumber();

  diagonalMatrix = MatrixNC::Zero(numPixels, numPixels);

  RBuffer2D buffer(cols, rows);
  buffer.assign(powerOf<2>(wavenum) +
    2 * wavenum / mParams.deltaZ() * mPhaseShifts->phaseShift(sliceIndex));

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      diagonalMatrix.diagonal()[y * cols + x] = buffer.pixel(x, y);

  if (verbosity >= 1)
  {
    std::cout << "   * A (" << sliceIndex << ") Time for calcDiagonalPart: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void PassLaplace::calcLaplaceMatrix(MatrixNC& laplaceMatrix)
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  if (verbosity >= 1)
    std::cout << "- calcLaplaceMatrix started." << std::endl;

  const int64_t cols = mParams.cols();
  const int64_t rows = mParams.rows();
  const int64_t numPixels = mParams.numPixels();

  const Real deltaKX = mParams.deltaKX();
  const Real deltaKY = mParams.deltaKY();

  // inverse Fourier trafo
  MatrixNC invFourierTrafoMatrix;
  calcFourierTrafoMatrix(mParams, invFourierTrafoMatrix, -1.0);
  if (verbosity >= 1)
    std::cout << "   * inverse Fourier trafo matrix calculated." << std::endl;

  DiagonalMatrixNC diagMatrix(numPixels);

  #pragma omp parallel for
  for (int64_t j = 0; j < rows; ++j)
  {
    // coordinates in reciprocal space
    const Real kY = deltaKY * (j > rows / 2 ? j - rows : j);
    const Real kYSqr = kY * kY;

    for (int64_t i = 0; i < cols; ++i)
    {
      const Real kX = deltaKX * (i > cols / 2 ? i - cols : i);
      const Real kXSqr = kX * kX;

      diagMatrix.diagonal()[j * cols + i] = (kXSqr + kYSqr) / numPixels;
    }
  }

  invFourierTrafoMatrix *= diagMatrix;

  if (verbosity >= 1)
    std::cout << "   * Diagonal matrix applied." << std::endl;

  MatrixNC fourierTrafoMatrix;
  calcFourierTrafoMatrix(mParams, fourierTrafoMatrix, +1.0);
  if (verbosity >= 1)
    std::cout << "   * Fourier trafo matrix calculated." << std::endl;
  laplaceMatrix.noalias() = invFourierTrafoMatrix * fourierTrafoMatrix;

  if (verbosity >= 1)
  {
    std::cout << "   * Fourier trafo matrix applied.\n"
                 " - Time for calcLaplaceMatrix: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void PassLaplace::run()
{
  using namespace std::chrono;
  auto startTime = high_resolution_clock::now();

  MatrixNC laplace;
  calcLaplaceMatrix(laplace);
  saveMatrix(laplace, matrixDir() + "/laplace.mat");

  if (verbosity >= 1)
  {
    std::cout << "Laplace matrix calculated: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

PassMatrixA2::PassMatrixA2(const MultisliceParameters& params,
                           const std::string& matrixDir,
                           const ConstSamplePtr& sample, size_t slice) :
  Pass(params, matrixDir), mMatrixA2(params, sample, matrixDir + "laplace.mat"),
  mSlice(slice)
{}

void PassMatrixA2::run()
{
  auto startTime = high_resolution_clock::now();

  VectorNR eigenvalues;
  MatrixNC eigenvectors;
  mMatrixA2.calc(mSlice, eigenvalues, eigenvectors);
  saveEVDecomposition(eigenvalues, eigenvectors, matrixName("A-", mSlice));

  if (verbosity >= 1)
  {
    std::cout << "EV Decomposition of Matrix2A calculated: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

Pass1::Pass1(const MultisliceParameters& params, const std::string& matrixDir,
             const SamplePtr& sample, size_t numParts, size_t part) :
  Pass(params, matrixDir), mMatrixA2(params, sample, matrixDir + "laplace.mat"),
  mPart(params, numParts, part)
{}

void Pass1::run()
{
  auto startTime = high_resolution_clock::now();

  if (verbosity >= 1)
    std::cout << "Pass1 started." << std::endl;

  for (size_t i = mPart.startSlice(); i < mPart.endSlice(); ++i)
  {
    VectorNR eigenvaluesSqr;
    MatrixNC eigenvectors;

    mMatrixA2.calc(i, eigenvaluesSqr, eigenvectors);
    saveEVDecomposition(eigenvaluesSqr, eigenvectors, matrixName("A-", i));

    if (verbosity >= 1)
    {
      std::cout << "   * (" << i << ") Runtime: "
                << duration(high_resolution_clock::now() - startTime)
                << std::endl;
    }
  }
}

} // namespace Aurora
