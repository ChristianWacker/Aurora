//--- Aurora/Clients/Matrix/Pass1.hpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_MATRIX_PASS1_HPP
#define AURORA_CLIENTS_MATRIX_PASS1_HPP

#include "AuroraLib/PhaseShift.hpp"
#include "Pass.hpp"

namespace Aurora
{

/// Class to claculate the matrix A^2 or more precisely its eigenvalue
/// decomposition.
class MatrixA2
{
public:
  MatrixA2(const MultisliceParameters& params, const ConstSamplePtr& sample,
           const std::string& filenameLaplaceMatrix);

  void calc(size_t slice, VectorNR& eigenvaluesSqr,
            MatrixNC& eigenvectors) const;

private:
  MultisliceParameters mParams;
  std::unique_ptr<PhaseShift> mPhaseShifts;
  std::string mFilenameLaplaceMatrix;

  /// Calculate the part of the matrix A2 that is diagonal.
  void calcDiagonalPart(Eigen::Ref<MatrixNC> diagonalMatrix,
                        size_t sliceIndex) const;
};

class PassLaplace : public Pass
{
public:
  using Pass::Pass;
  void run() override;

private:
  void calcLaplaceMatrix(MatrixNC& laplaceMatrix);
};

class PassMatrixA2 : public Pass
{
public:
  PassMatrixA2(const MultisliceParameters& params, const std::string& matrixDir,
               const ConstSamplePtr& sample, size_t slice);

  void run() override;

private:
  MatrixA2 mMatrixA2;
  size_t mSlice;
};

/// This is the first pass. It calculates the eigenvalue decomposition for all
/// matrices A^2.
class Pass1 : public Pass
{
public:
  Pass1(const MultisliceParameters& params, const std::string& matrixDir,
        const SamplePtr& sample, size_t numParts, size_t part);

  void run() override;

private:
  MatrixA2 mMatrixA2;
  Part mPart;
};

} // namespace Aurora

#endif
