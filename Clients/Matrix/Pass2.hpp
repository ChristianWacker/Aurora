//--- Aurora/Clients/Matrix/Pass2.hpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_MATRIX_PASS2_HPP
#define AURORA_CLIENTS_MATRIX_PASS2_HPP

#include "Pass.hpp"
#include "Pass1.hpp"

namespace Aurora
{

class Pass2 : public Pass
{
public:
  Pass2(const MultisliceParameters& params, const std::string& matrixDir,
        const ConstSamplePtr& sample, size_t numParts, size_t part) :
    Pass(params, matrixDir), mPart(params, numParts, part)
  {
    if (sample)
    {
      mMatrixA2 = std::make_unique<MatrixA2>(params, sample,
                                             matrixDir + "laplace.mat");
    }
  }

protected:
  Part mPart;
  std::unique_ptr<MatrixA2> mMatrixA2;
};

/// Combine the matrices A for the Matrix-Full algorithm.
class Pass2F : public Pass2
{
public:
  using Pass2::Pass2;
  void run() override;
};

class Pass2FLarge : public Pass2F
{
public:
  using Pass2F::Pass2F;
  void run() override;
};

/// Combine the matrices A for the matrix algorithm
class Pass2T : public Pass2
{
public:
  using Pass2::Pass2;
  void run() override;
};

class Pass2TLarge : public Pass2T
{
public:
  using Pass2T::Pass2T;
  void run() override;
};

/// Combine the matrices A for the matrix algorithm without bandlimit
class Pass2TSLarge : public Pass2
{
public:
  using Pass2::Pass2;
  void run() override;
};

/// Combine the matrices A for the matrix algorithm with bandlimit
class Pass2TSBandlimitLarge : public Pass2TSLarge
{
public:
  using Pass2TSLarge::Pass2TSLarge;
  void run() override;
};

} // namespace Aurora

#endif
