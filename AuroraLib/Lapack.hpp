//--- Aurora/AuroraLib/Lapack.hpp ----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_LAPACK_HPP
#define AURORA_AURORA_LIB_LAPACK_HPP

#include "AuroraLib/Prerequisites.hpp"

#include "AuroraLib/Complex.hpp"
#include "AuroraLib/Exception.hpp"
#include "AuroraLib/MatrixSupport.hpp"

#include <iostream>

#if (AURORA_USE_LAPACK == 1)

namespace Aurora
{

namespace Lapack
{

typedef MKL_INT Int;
typedef MKL_Complex16 Complex;

/// Calculate the eigenvalues and eigenvectors using Lapack. We try to create an
/// interface is reasonably compatible with the Eigen library.
class SelfAdjointEigenSolver
{
public:
  SelfAdjointEigenSolver(const MatrixNC& mat) :
    mInfo(Eigen::ComputationInfo::NoConvergence)
  {
    compute(mat);
  }

  Eigen::ComputationInfo info() const { return mInfo; }
  const VectorNR& eigenvalues() const { return mEigenvalues; }
  const MatrixNC& eigenvectors() const { return mEigenvectors; }

private:
  void compute(const MatrixNC& mat)
  {
    // Although the function zheev has a much simpler interface the Intel
    // documentation claims that zheevr is faster and uses less workspace.

    mMat = mat;

    assert(mMat.rows() == mMat.cols() && "Matrix must be square.");

    // We emulate the Eigen interface
    mInfo = Eigen::ComputationInfo::NoConvergence;

    const Int dimension = mMat.rows();
    double epsilon = 1e-6;
    Int numEigenvaluesFound = 0;

    // allocate space for the eigenvalues and eigenvectors
    mEigenvalues.resize(dimension);
    mEigenvectors.resize(dimension, dimension);

    // unused but needed dummy variables
    double dummyD = 0.0;
    Int dummyI = 0;

    auto iSupport = std::make_unique<Int[]>(2 * dimension);

    // prepare everything for a workspace query
    auto cWork = std::make_unique<Complex[]>(1);
    auto rWork = std::make_unique<double[]>(1);
    auto iWork = std::make_unique<Int[]>(1);

    Int lengthCWork = -1;
    Int lengthRWork = -1;
    Int lengthIWork = -1;
    Int info = 0;

    // issue the workspace query
    zheevr_("V", "A", "U", &dimension,
      reinterpret_cast<Complex*>(mMat.data()), &dimension, &dummyD, &dummyD,
      &dummyI, &dummyI, &epsilon, &numEigenvaluesFound, mEigenvalues.data(),
      reinterpret_cast<Complex*>(mEigenvectors.data()), &dimension,
      iSupport.get(), cWork.get(), &lengthCWork, rWork.get(), &lengthRWork,
      iWork.get(), &lengthIWork, &info);

    if (0 != info)
      AURORA_THROW(EUnknown, "zheevr failure (query): " + std::to_string(info));

    // extract the workspace size
    lengthCWork = static_cast<Int>(cWork[0].real);
    lengthRWork = static_cast<Int>(rWork[0]);
    lengthIWork = iWork[0];

    // allocate the workspace
    cWork = std::make_unique<Complex[]>(lengthCWork);
    rWork = std::make_unique<double[]>(lengthRWork);
    iWork = std::make_unique<Int[]>(lengthIWork);

    // execute the actual calculation
    zheevr_("V", "A", "L", &dimension,
      reinterpret_cast<Complex*>(mMat.data()), &dimension, &dummyD, &dummyD,
      &dummyI, &dummyI, &epsilon, &numEigenvaluesFound, mEigenvalues.data(),
      reinterpret_cast<Complex*>(mEigenvectors.data()), &dimension,
      iSupport.get(), cWork.get(), &lengthCWork, rWork.get(), &lengthRWork,
      iWork.get(), &lengthIWork, &info);

    if ((0 != info) || (dimension != numEigenvaluesFound))
    {
      AURORA_THROW(EUnknown,
        "zheevr failure (calc): " + std::to_string(info) +
        ", numEigenvaluesFound: " + std::to_string(numEigenvaluesFound));
    }

    mInfo = Eigen::ComputationInfo::Success;
  }

  Eigen::ComputationInfo mInfo;

  VectorNR mEigenvalues;
  MatrixNC mEigenvectors;
  // We store a copy of the matrix and not a reference as the matrix wil be
  // overwritten by LAPACK
  MatrixNC mMat;
};

class PartialPivLU
{
public:
  typedef Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, Int>
    PermutationMatrix;

  PartialPivLU() = default;

  PartialPivLU(const MatrixNC& matrix) :
    mLU(matrix.rows(), matrix.rows()),
    mPermutation(matrix.rows())
  {
    compute(matrix);
  }

  PartialPivLU& compute(const MatrixNC& matrix)
  {
    mLU = matrix;

    const Int dimension = mLU.rows();

    mPermutation.resize(dimension);

    Int info = 0;

    zgetrf_(&dimension, &dimension,
      reinterpret_cast<Complex*>(mLU.data()), &dimension,
      mPermutation.indices().data(), &info);

    if (0 != info)
      AURORA_THROW(EUnknown, "zgetrf failure: " + std::to_string(info));

    return *this;
  }

  inline const PermutationMatrix& permutationP() const
  {
    return mPermutation;
  }

  inline const MatrixNC& matrixLU() const
  {
    return mLU;
  }

  VectorNC solve(const VectorNC& b)
  {
    const Int dimension = mLU.rows();
    // number of columns of b
    const Int nhrs = 1;

    VectorNC x = b;

    Int info = 0;

    zgetrs("N", &dimension, &nhrs,
      reinterpret_cast<Complex*>(mLU.data()), &dimension,
      mPermutation.indices().data(),
      reinterpret_cast<Complex*>(x.data()), &dimension, &info);

    if (0 != info)
      AURORA_THROW(EUnknown, "zgetrs failure: " + std::to_string(info));

    return x;
  }

private:
  MatrixNC mLU;
  PermutationMatrix mPermutation;
};

} // namespace Lapack

typedef Lapack::SelfAdjointEigenSolver SelfAdjointEigenSolver;
typedef Lapack::PartialPivLU PartialPivLU;

} // namespace Aurora

#else

namespace Aurora
{
  // Fallback to Eigen
  typedef Eigen::SelfAdjointEigenSolver<MatrixNC> SelfAdjointEigenSolver;
  typedef Eigen::PartialPivLU<MatrixNC> PartialPivLU;
} // namespace Aurora

#endif // (AURORA_USE_LAPACK == 1)

#endif
