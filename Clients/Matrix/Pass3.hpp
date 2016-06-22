//--- Aurora/Clients/Matrix/Pass3.hpp ------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_MATRIX_PASS3_HPP
#define AURORA_CLIENTS_MATRIX_PASS3_HPP

#include "Pass.hpp"

namespace Aurora
{

/// Class for the final combination of the transfer matrices.
/// The singles results from the previous steps will be combined by this class.
class Pass3 : public Pass
{
public:
  /// @param matrixPrefix should be "F", "T" or "S"
  Pass3(const MultisliceParameters& params, const std::string& matrixDir,
        size_t numParts, const std::string& matrixPrefix, Direction direction) :
    Pass(params, matrixDir), mNumParts(numParts), mMatrixPrefix(matrixPrefix),
    mDirection(direction)
  {
    assert(mNumParts >= 1 && "Invalid number of parts.");
  }

  void run() override;

protected:
  size_t mNumParts;
  std::string mMatrixPrefix;
  Direction mDirection;
};

class Pass3Large : public Pass3
{
public:
  using Pass3::Pass3;
  void run() override;
};

} // namespace Aurora

#endif
