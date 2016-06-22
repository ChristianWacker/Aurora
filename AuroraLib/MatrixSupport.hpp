//--- Aurora/AuroraLib/MatrixSupport.cpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MATRIX_SUPPORT_HPP
#define AURORA_AURORA_LIB_MATRIX_SUPPORT_HPP

#include "AuroraLib/Complex.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace Aurora
{

typedef Eigen::Matrix<Real, 3, 3> Matrix3x3R;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> MatrixNC;
typedef Eigen::DiagonalMatrix<Real, Eigen::Dynamic> DiagonalMatrixNR;
typedef Eigen::DiagonalMatrix<Complex, Eigen::Dynamic> DiagonalMatrixNC;

typedef Eigen::Matrix<Real, 3, 1> Vector3R;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorNR;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> VectorNC;

AURORA_API void loadMatrix(const std::string& filename, MatrixNC& matrix);
AURORA_API void saveMatrix(const MatrixNC& mat, const std::string& filename);

AURORA_API void loadEVDecomposition(const std::string& filename,
                                    VectorNR& eigenvalues,
                                    MatrixNC& eigenvectors);
AURORA_API void loadEigenvectors(const std::string& filename,
                                 MatrixNC& eigenvectors);
AURORA_API void loadEigenvalues(const std::string& filename,
                                VectorNR& eigenvalues);
AURORA_API void saveEVDecomposition(const VectorNR& eigenvalues,
                                    const MatrixNC& eigenvectors,
                                    const std::string& filename);

AURORA_API void inspect(const std::string& matrixName, int submatrix = 0);

AURORA_API void sqrtEigenValues(const VectorNR& eigenValuesSqr,
                                VectorNC& eigenValues);

/// Converts a Vector into CBuffer2D. The buffer must have the correct
/// dimensions.
AURORA_API void vectorToBuffer(const VectorNC& vec, CBuffer2D& buffer);
AURORA_API void bufferToVector(const CBuffer2D& buffer, VectorNC& vec);


AURORA_API void calcFourierTrafoMatrix(const MultisliceParameters& params,
                                       MatrixNC& fourierTrafoMatrix, Real sign);

AURORA_API void multiplyLarge(const std::string& src1, const std::string& src2,
                              const std::string& dst);

AURORA_API void removeMatricesLarge(const std::string& baseName);

} // namespace Aurora

#endif
