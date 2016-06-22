//--- Aurora/AuroraLib/SmallVector.hpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_SMALL_VECTOR_HPP
#define AURORA_AURORA_LIB_SMALL_VECTOR_HPP

#include "AuroraLib/Prerequisites.hpp"

namespace Aurora
{

/// @brief
///  A vector class optimized for a small number of elements.
template<class T, size_t n>
class SmallVector
{
public:
  SmallVector(size_t size) :
    mSize(size)
  {
    if (mSize > n)
      mData = new T[mSize];
    else
      mData = mStaticData;
  }

  ~SmallVector()
  {
    if (mSize > n)
      delete[] mData;
  }

  T* get() const
  {
    return mData;
  }

  T& operator[](size_t index)
  {
    return mData[index];
  }

  size_t size() const
  {
    return mSize;
  }

private:
  T mStaticData[n];
  T* mData;
  size_t mSize;
};

} // namespace Aurora

#endif
