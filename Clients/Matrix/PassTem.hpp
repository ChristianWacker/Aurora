//--- Aurora/Clients/Matrix/PassTem.hpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_MATRIX_PASS_TEM_HPP
#define AURORA_CLIENTS_MATRIX_PASS_TEM_HPP

#include "Pass.hpp"

namespace Aurora
{

class PassTem : public PassFinal
{
public:
  using PassFinal::PassFinal;

  void run() override;
};

class PassTemLarge : public PassTem
{
public:
  using PassTem::PassTem;

  void run() override;
};

class PassTemIterative : public PassTem
{
public:
  PassTemIterative(const MultisliceParameters& params,
    const std::string& matrixDir, size_t power, const std::string& outTemplate,
    size_t minIterations, size_t maxIterations);

  void run() override;

private:
  size_t mMinIterations, mMaxIterations;
};

} // namespace Aurora

#endif
