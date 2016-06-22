//--- Aurora/Clients/Matrix/Pass.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_MATRIX_PASS_HPP
#define AURORA_CLIENTS_MATRIX_PASS_HPP

#include "AuroraLib/MatrixSupport.hpp"
#include "AuroraLib/MultisliceParameters.hpp"

namespace Aurora
{

/// Class representing a single processing step for the calculation of the
/// final transfer matrix.
class Pass
{
public:
  Pass(const MultisliceParameters& params, const std::string& matrixDir) :
    mParams(params), mMatrixDir(matrixDir)
  { }

  virtual ~Pass() { }
  virtual void run() = 0;

  std::string matrixDir() const { return mMatrixDir; }

protected:
    MultisliceParameters mParams;

  std::string matrixName(const std::string& prefix,
                       size_t index = std::numeric_limits<size_t>::max()) const;
  std::string tempName(const std::string& prefix, size_t index) const;

  void dumpHistogram(const MatrixNC& matrix, const std::string& filename);

private:
  const std::string mMatrixDir;
};

/// Base class for all passes that output waves
class PassWaveOut
{
public:
  PassWaveOut(const std::string& outTemplate) :
    mOutTemplate(outTemplate)
  { }

  void saveWave(const CBuffer2D& wave, const Buffer2DInfo& params,
                const std::string& infix, size_t index) const;

private:
  const std::string mOutTemplate;
};

class PassFinal : public Pass, public PassWaveOut
{
public:
  PassFinal(const MultisliceParameters& params, const std::string& matrixDir,
    size_t power, const std::string& outTemplate);

protected:
  size_t mPower;
};

class Part
{
public:
  Part(const MultisliceParameters& params, size_t numParts, size_t part);

  size_t numParts() const { return mNumParts; }
  operator size_t() const { return mPart; }
  size_t startSlice() const { return mStartSlice; }
  size_t endSlice() const { return mEndSlice; }

private:
  const size_t mNumParts;
  const size_t mPart;

  size_t mStartSlice;
  size_t mEndSlice;
};

} // namespace Aurora

#endif
