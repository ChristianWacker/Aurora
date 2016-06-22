//--- Aurora/AuroraLib/Histogram.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_HISTOGRAM_HPP
#define AURORA_AURORA_LIB_HISTOGRAM_HPP

#include "AuroraLib/Math.hpp"
#include "AuroraLib/Prerequisites.hpp"
#include "AuroraLib/Utils.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace Aurora
{

/// This class represents a one-dimensional histogram
class AURORA_API Histogram
{
public:
  Histogram(Real minValue, Real maxValue, size_t numBins) :
    mMinValue(minValue), mMaxValue(maxValue), mBins(numBins, 0)
  {
    assert(minValue != maxValue);
  }

  // Returns the number of bins
  size_t numBins() const
  {
    return mBins.size();
  }

  size_t operator[](size_t bin) const 
  {
    return mBins[bin];
  }

  void insert(Real value)
  {
    ++mBins[valueToBin(value)];
  }

  void reset()
  {
    std::fill(mBins.begin(), mBins.end(), 0);
  }

  Real binToValue(size_t bin) const
  {
    return static_cast<Real>(bin) / static_cast<Real>(mBins.size()) *
      (mMaxValue - mMinValue) + mMinValue
      + R(0.5) * (mMaxValue - mMinValue) / static_cast<Real>(mBins.size());
  }

  size_t valueToBin(Real value) const
  {
    Real scaledValue = (value - mMinValue) / (mMaxValue - mMinValue);
    return clamp<int64_t>(static_cast<int64_t>(scaledValue * mBins.size()),
      0, mBins.size() - 1);
  }

  Real minValue() const
  {
    return mMinValue;
  }

  Real maxValue() const
  {
    return mMaxValue;
  }

  /// Returns the bin with the greatest number of entries
  size_t mode() const
  {
    return std::max_element(mBins.begin(), mBins.end()) - mBins.begin();
  }

  size_t numEntries() const
  {
    return std::accumulate(mBins.begin(), mBins.end(), 0);
  }

  Real mean() const
  {
    Real sum = 0;
    for (size_t bin = 0, numBins = mBins.size(); bin < numBins; ++bin)
      sum += binToValue(bin) * mBins[bin];

    return sum / numEntries();
  }

  Real variance() const
  {
    Real sum = 0;
    const size_t numBins = mBins.size();

    const Real mean = this->mean();
    for (size_t bin = 0; bin < numBins; ++bin)
      sum += powerOf<2>(binToValue(bin) - mean) * mBins[bin];

    return sum / numEntries();
  }

  void saveAsText(const std::string& filename,
                  Real trueMean = nan(""), Real trueStddev = nan("")) const;

private:
  const Real mMinValue;
  const Real mMaxValue;
  std::vector<size_t> mBins;
};

} // namespace Aurora

#endif
