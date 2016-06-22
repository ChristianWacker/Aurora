//--- Aurora/AuroraLib/MatrixSupport.cpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/MatrixSupport.hpp"

#include "AuroraLib/BinaryIStream.hpp"
#include "AuroraLib/BinaryOStream.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/MultisliceParameters.hpp"


#include <chrono>
#include <iostream>

namespace Aurora
{

using namespace std::chrono;

namespace
{
  struct MatrixFileHeader
  {
    uint32_t version = 0;
    uint32_t type = 0;
    int64_t  cols = 0;
    int64_t  rows = 0;
  };

  const uint32_t matrixFileVersion = 2;

  const uint32_t matrixTypeComplex = 1;
  const uint32_t matrixTypeSelfAdjointDecomposition = 2;
}

void loadMatrix(const std::string& filename, MatrixNC& matrix)
{
  auto startTime = high_resolution_clock::now();
  if (verbosity >= 1)
    std::cout << "Loading matrix \"" << filename << "\"... " << std::flush;

  BinaryIStream stream = BinaryIStream::openFile(filename);

  MatrixFileHeader header;
  stream.read<uint32_t>(&header.version, 2);
  stream.read<int64_t>(&header.cols, 2);

  if ((matrixFileVersion != header.version) ||
      (matrixTypeComplex != header.type) ||
      (1 > header.cols) || (1 > header.rows))
    AURORA_THROW(EInOut, "Cannot load \"" + filename + "\".");

  const size_t numElements = header.cols * header.rows;
  matrix.resize(header.cols, header.rows);
  stream.read<Complex>(matrix.data(), numElements);

  if (verbosity >= 1)
  {
    std::cout << "done. Time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void saveMatrix(const MatrixNC& matrix, const std::string& filename)
{
  auto startTime = high_resolution_clock::now();
  if (verbosity >= 1)
    std::cout << "Saving matrix \"" << filename << "\"... " << std::flush;

  BinaryOStream stream = BinaryOStream::createFile(filename);
  stream.write<uint32_t>(matrixFileVersion);
  stream.write<uint32_t>(matrixTypeComplex);
  stream.write<int64_t>(matrix.cols());
  stream.write<int64_t>(matrix.rows());

  stream.write<Complex>(matrix.data(), matrix.cols() * matrix.rows());

  if (verbosity >= 1)
  {
    std::cout << "done. Time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void loadEVDecomposition(const std::string& filename, VectorNR& eigenvalues,
                         MatrixNC& eigenvectors)
{
  auto startTime = high_resolution_clock::now();
  if (verbosity >= 1)
  {
    std::cout << "Loading decomposition \"" << filename << "\"... "
              << std::flush;
  }

  BinaryIStream stream = BinaryIStream::openFile(filename);

  MatrixFileHeader header;
  stream.read<uint32_t>(&header.version, 2);
  stream.read<int64_t>(&header.cols, 2);

  if ((matrixFileVersion != header.version) ||
      (matrixTypeSelfAdjointDecomposition != header.type) ||
      (1 > header.cols) || (1 > header.rows))
    AURORA_THROW(EInOut, "Cannot load \"" + filename + "\".");

  eigenvalues.resize(header.rows);
  stream.read<Real>(eigenvalues.data(), header.rows);

  const size_t numElements = header.cols * header.rows;
  eigenvectors.resize(header.rows, header.cols);
  stream.read<Complex>(eigenvectors.data(), numElements);

  if (verbosity >= 1)
  {
    std::cout << "done. Time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void loadEigenvectors(const std::string& filename, MatrixNC& eigenvectors)
{
  VectorNR eigenvalues;
  loadEVDecomposition(filename, eigenvalues, eigenvectors);
}

void loadEigenvalues(const std::string& filename, VectorNR& eigenvalues)
{
  auto startTime = high_resolution_clock::now();
  if (verbosity >= 1)
  {
    std::cout << "Loading eigenvalues \"" << filename << "\"... "
              << std::flush;
  }

  BinaryIStream stream = BinaryIStream::openFile(filename);

  MatrixFileHeader header;
  stream.read<uint32_t>(&header.version, 2);
  stream.read<int64_t>(&header.cols, 2);

  if ((matrixFileVersion != header.version) ||
      (matrixTypeSelfAdjointDecomposition != header.type) ||
      (1 > header.cols) || (1 > header.rows))
    AURORA_THROW(EInOut, "Cannot load \"" + filename + "\".");

  eigenvalues.resize(header.rows);
  stream.read<Real>(eigenvalues.data(), header.rows);

  if (verbosity >= 1)
  {
    std::cout << "done. Time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void saveEVDecomposition(const VectorNR& eigenvalues,
                       const MatrixNC& eigenvectors,
                       const std::string& filename)
{
  auto startTime = high_resolution_clock::now();
  if (verbosity >= 1)
  {
    std::cout << "Saving decomposition \"" << filename << "\"... "
              << std::flush;
  }

  assert(eigenvalues.rows() == eigenvectors.rows());

  BinaryOStream stream = BinaryOStream::createFile(filename);
  stream.write<uint32_t>(matrixFileVersion);
  stream.write<uint32_t>(matrixTypeSelfAdjointDecomposition);
  stream.write<int64_t>(eigenvectors.cols());
  stream.write<int64_t>(eigenvectors.rows());

  stream.write<Real>(eigenvalues.data(), eigenvectors.rows());
  stream.write<Complex>(eigenvectors.data(),
                        eigenvectors.rows() * eigenvectors.cols());

  if (verbosity >= 1)
  {
    std::cout << "done. Time needed: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void inspect(const std::string& matrixName, int submatrixIndex)
{
  BinaryIStream stream = BinaryIStream::openFile(matrixName);

  MatrixFileHeader header;
  stream.read<uint32_t>(&header.version, 2);
  stream.read<int64_t>(&header.cols, 2);

  if (matrixFileVersion != header.version)
    AURORA_THROW(EInOut, "Cannot load \"" + matrixName + "\".");

  VectorNC eigenValues;

  if (matrixTypeSelfAdjointDecomposition == header.type)
  {
    VectorNR eigenValuesFile;
    loadEigenvalues(matrixName, eigenValuesFile);

    eigenValues = eigenValuesFile.cast<Complex>();
  }
  else if (matrixTypeComplex == header.type)
  {
    MatrixNC matrix;
    loadMatrix(matrixName, matrix);

    const int64_t numPixels = matrix.cols() / 2;

    MatrixNC submatrix;
    switch (submatrixIndex)
    {
    case 1:
      submatrix = matrix.block(0, 0, numPixels, numPixels);
      break;
    case 2:
      submatrix = matrix.block(0, numPixels, numPixels, numPixels);
      break;
    case 3:
      submatrix = matrix.block(numPixels, 0, numPixels, numPixels);
      break;
    case 4:
      submatrix = matrix.block(numPixels, numPixels, numPixels, numPixels);
      break;
    default:
      submatrix = matrix;
    }

    // the matrix t11
    Eigen::ComplexEigenSolver<MatrixNC> decomposition(submatrix, false);
    if (decomposition.info() != Eigen::ComputationInfo::Success)
      AURORA_THROW(EUnknown, "Decomposition could not be calculated.");

    eigenValues = decomposition.eigenvalues();
  }
  else
  {
    AURORA_THROW(EInOut, "Cannot load \"" + matrixName + "\".");
  }

  int64_t numPixels = eigenValues.rows();
  std::cerr << "numPixels: " << numPixels << '\n';

  std::cout << "# Filename: " << matrixName << "\n"
               "# Real             Imag             Norm\n";
  std::cout.precision(12);
  for (int64_t i = 0; i < numPixels; ++i)
  {
    Complex c = eigenValues[i];
    std::cout << c.real() << ' ' << c.imag() << ' ' << std::norm(c) << '\n';
  }
}

void sqrtEigenValues(const VectorNR& eigenValuesSqr, VectorNC& eigenValues)
{
  const int64_t numEigenValues = eigenValuesSqr.rows();
  eigenValues.resize(numEigenValues);

  #pragma omp parallel for
  for (int64_t i = 0; i < numEigenValues; ++i)
  {
    if (eigenValuesSqr[i] > R(0.0))
      eigenValues[i] = std::sqrt(eigenValuesSqr[i]);
    else
      eigenValues[i] = Complex(R(0.0), std::sqrt(-eigenValuesSqr[i]));
  }
}

void vectorToBuffer(const VectorNC& vec, CBuffer2D& buffer)
{
  const int64_t cols = buffer.cols();
  const int64_t rows = buffer.rows();

  assert(vec.rows() == cols * rows);

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      buffer.pixel(x, y) = vec[y * cols + x];
}

void bufferToVector(const CBuffer2D& buffer, VectorNC& vec)
{
  const int64_t cols = buffer.cols();
  const int64_t rows = buffer.rows();

  vec.resize(cols * rows);

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      vec[y * cols + x] = buffer.pixel(x, y);
}

void calcFourierTrafoMatrix(const MultisliceParameters& params,
                            MatrixNC& fourierTrafoMatrix, Real sign)
{
  auto startTime = high_resolution_clock::now();

  const int64_t cols = params.cols();
  const int64_t rows = params.rows();
  const int64_t numPixels = cols * rows;

  fourierTrafoMatrix.resize(numPixels, numPixels);

  #pragma omp parallel for
  for (int64_t outY = 0; outY < rows; ++outY)
  {
    for (int64_t outX = 0; outX < cols; ++outX)
    {
      const int64_t outIdx = outY * cols + outX;
      for (int64_t inY = 0; inY < rows; ++inY)
      {
        for (int64_t inX = 0; inX < cols; ++inX)
        {
          const int64_t inIdx = inY * cols + inX;

          Real phase = Math::twoPi * (outX * inX / static_cast<Real>(cols) +
                                      outY * inY / static_cast<Real>(rows));
          fourierTrafoMatrix(outIdx, inIdx) = std::polar(R(1.0), sign * phase);
        }
      }
    }
  }

  if (verbosity >= 1)
  {
    std::cout << "   * Time for calcFourierTrafoMatrix: "
              << duration(high_resolution_clock::now() - startTime)
              << std::endl;
  }
}

void multiplyLarge(const std::string& src1, const std::string& src2,
                   const std::string& dst)
{
  {
    MatrixNC c11;
    {
      MatrixNC a11;
      loadMatrix(src1 + "-11.mat", a11);
      MatrixNC b11;
      loadMatrix(src2 + "-11.mat", b11);
      c11.noalias() = a11 * b11;
    }
    {
      MatrixNC a12;
      loadMatrix(src1 + "-12.mat", a12);
      MatrixNC b21;
      loadMatrix(src2 + "-21.mat", b21);
      c11.noalias() += a12 * b21;
    }
    saveMatrix(c11, dst + "-11.mat");
  }

  {
    MatrixNC c12;
    {
      MatrixNC a11;
      loadMatrix(src1 + "-11.mat", a11);
      MatrixNC b12;
      loadMatrix(src2 + "-12.mat", b12);
      c12.noalias() = a11 * b12;
    }
    {
      MatrixNC a12;
      loadMatrix(src1 + "-12.mat", a12);
      MatrixNC b22;
      loadMatrix(src2 + "-22.mat", b22);
      c12.noalias() += a12 * b22;
    }
    saveMatrix(c12, dst + "-12.mat");
  }

  {
    MatrixNC c21;
    {
      MatrixNC a21;
      loadMatrix(src1 + "-21.mat", a21);
      MatrixNC b11;
      loadMatrix(src2 + "-11.mat", b11);
      c21.noalias() = a21 * b11;
    }
    {
      MatrixNC a22;
      loadMatrix(src1 + "-22.mat", a22);
      MatrixNC b21;
      loadMatrix(src2 + "-21.mat", b21);
      c21.noalias() += a22 * b21;
    }
    saveMatrix(c21, dst + "-21.mat");
  }

  {
    MatrixNC c22;
    {
      MatrixNC a21;
      loadMatrix(src1 + "-21.mat", a21);
      MatrixNC b12;
      loadMatrix(src2 + "-12.mat", b12);
      c22.noalias() = a21 * b12;
    }
    {
      MatrixNC a22;
      loadMatrix(src1 + "-22.mat", a22);
      MatrixNC b22;
      loadMatrix(src2 + "-22.mat", b22);
      c22.noalias() += a22 * b22;
    }
    saveMatrix(c22, dst + "-22.mat");
  }
}

void removeMatricesLarge(const std::string& baseName)
{
  remove((baseName + "-11.mat").c_str());
  remove((baseName + "-12.mat").c_str());
  remove((baseName + "-21.mat").c_str());
  remove((baseName + "-22.mat").c_str());
}

} // namespace Aurora
