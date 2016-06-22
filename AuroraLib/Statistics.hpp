//--- Aurora/AuroraLib/Statistics.hpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_STATISTICS_HPP
#define AURORA_AURORA_LIB_STATISTICS_HPP

#include "AuroraLib/Math.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace Aurora
{

template<class T>
class Statistics
{
public:
  T mean() const
  {
    return accumulate(mData.begin(), mData.end(), T(0)) / mData.size();
  }

  T median() const
  {
    return mData[mData.size() / 2];
  }

  T variance() const
  {
    T mu(mean());
    auto f = [mu] (T sum, T value) { return sum + powerOf<2>(value - mu); };
    return accumulate(mData.begin(), mData.end(), T(0), f) / (mData.size() - 1);
  }

  T standardError() const
  {
    return sqrt(variance() / mData.size());
  }

  T min() const
  {
    return *std::min_element(mData.begin(), mData.end());
  }

  T max() const
  {
    return *std::max_element(mData.begin(), mData.end());
  }

  void add(T value)
  {
    mData.push_back(value);
  }

private:
  std::vector<T> mData;
};

} // namespace Aurora

#endif
