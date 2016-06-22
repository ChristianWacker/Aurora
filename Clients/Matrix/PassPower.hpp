//--- Aurora/Clients/Matrix/PassPower.hpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_MATRIX_PASS_POWER_HPP
#define AURORA_CLIENTS_MATRIX_PASS_POWER_HPP

#include "Pass.hpp"

namespace Aurora
{

/// Powers of the matrix will be calculated for the simulation of thicker
/// samples. The final results is then sent to the processMatrix() function that
/// needs to be specialized in the derived classes.
class PassPower : public Pass
{
public:
  PassPower(const MultisliceParameters& params, const std::string& matrixDir,
    size_t power, const std::string& matrixPrefix, size_t startIndex = 1,
    size_t increment = 1);

  void run() override;

protected:
  size_t mPower;
  std::string mMatrixPrefix;
  size_t mStartIndex;
  size_t mIncrement;

  virtual void processMatrix(size_t power, const MatrixNC& mat) = 0;
};

/// Class to calculate the exit wave from forward scattering only (Matrix-Full).
class PassPowerF : public PassFinal
{
public:
  using PassFinal::PassFinal;

private:
  virtual void run() override;
};

/// Matrix-Algorithm. This class calculates and saves the three matrices:
/// * T_11
/// * T_12
/// * T_22^-1 * T_21
class PassPowerT : public PassPower
{
public:
  PassPowerT(const MultisliceParameters& params, const std::string& matrixDir,
             size_t power);

private:
  void processMatrix(size_t power, const MatrixNC& mat) override;
};

class PassPowerLarge : public PassPower
{
public:
  using PassPower::PassPower;

private:
  void run() override;
  void processMatrix(size_t, const MatrixNC&) final { }
};

} // namespace Aurora

#endif
