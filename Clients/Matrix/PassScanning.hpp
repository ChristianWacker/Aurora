//--- Aurora/Clients/Matrix/PassScanning.hpp -----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_MATRIX_PASS_SCANNING_HPP
#define AURORA_CLIENTS_MATRIX_PASS_SCANNING_HPP

#include "Matrix.hpp"

#include "AuroraLib/MultisliceParameters.hpp"

#include "Pass.hpp"

namespace Aurora
{

struct Signal
{
  enum Type
  {
    none    = 0,
    probe   = 1 << 0,
    psiBack = 1 << 1,
    psiEx   = 1 << 2,
    sem     = 1 << 3,
    stem    = 1 << 4,
    all     = 0xffffffff
  };
};

class PassScanning : public PassFinal
{
public:
  PassScanning(const MultisliceParameters& params, const std::string& matrixDir,
               size_t power, const std::string& outTemplate,
               const std::shared_ptr<StemDetector>& stemDetector);

  void run() override;

  void setSignals(Signal::Type signal)
  {
    mSignals |= signal;
  }

  void unsetSignals(Signal::Type signal)
  {
    mSignals &= ~signal;
  }

  bool signal(Signal::Type signal)
  {
    return 0 != (mSignals & signal);
  }

protected:
  uint32_t mSignals = Signal::none;

  std::shared_ptr<StemDetector> mStemDetector;
};

class PassScanningIterative : public PassScanning
{
public:
  PassScanningIterative(const MultisliceParameters& params,
    const std::string& matrixDir, size_t power, const std::string& outTemplate,
    const std::shared_ptr<StemDetector>& stemDetector, int64_t scanCols,
    int64_t scanRows, size_t minIterations, size_t maxIterations);

  void run() override;

private:
  int64_t mScanCols, mScanRows;
  size_t mMinIterations, mMaxIterations;
};

} // namespace Aurora

#endif
