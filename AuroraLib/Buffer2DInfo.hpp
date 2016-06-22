//--- Aurora/AuroraLib/Buffer2DInfo.hpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_BUFFER_2D_INFO_HPP
#define AURORA_AURORA_LIB_BUFFER_2D_INFO_HPP

#include "AuroraLib/Prerequisites.hpp"

#include "AuroraLib/Energy.hpp"

namespace Aurora
{

class AURORA_API Buffer2DInfoBase
{
public:
  Buffer2DInfoBase() = default;

  Buffer2DInfoBase(int64_t cols, int64_t rows, Real width, Real height) :
    mCols(cols), mRows(rows), mWidth(width), mHeight(height)
  {
    update();
  }

  /// Returns the aspect ratio of simulation box in x and y direction
  Real aspectRatio() const { return width() / height(); }

  /// Returns the number of pixel columns in the buffer.
  int64_t cols() const { return mCols; }

  /// Returns the distance of two pixels in Fourier space in Angstrom^-1
  Real deltaKX() const { return mDeltaKX; }

  /// Returns the distance of two pixel in Fourier space in Angstrom^-1
  Real deltaKY() const { return mDeltaKY; }

  /// Returns the distance of two pixels in x direction in Angstrom
  Real deltaX() const { return mDeltaX; }

  /// Returns the distance of two pixels in y direction in Angstrom
  Real deltaY() const { return mDeltaY; }

  /// Return the dimension of the simulation box in y direction in Angstrom
  Real height() const { return mHeight; }

  Real kXNyquist() const { return mKXNyquist; }

  Real kYNyquist() const { return mKYNyquist; }

  template<class T = int64_t>
  T numPixels() const { return static_cast<T>(mCols * mRows); }

  Real pixelAspectRatio() const { return deltaX() / deltaY(); }

  /// Returns the number of pixel rows in the buffer
  int64_t rows() const { return mRows; }

  /// Returns the dimension of the simulation box in x direction in Angstrom
  Real width() const { return mWidth; }

  bool operator==(const Buffer2DInfoBase& other) const
  {
    return (mCols == other.mCols) && (mRows == other.mRows) &&
           (mWidth == other.mWidth) && (mHeight == other.mHeight);
  }

  bool operator!=(const Buffer2DInfoBase& other) const
  {
    return !operator==(other);
  }

  void setResolution(int64_t cols, int64_t rows)
  {
    mRows = rows;
    mCols = cols;
    update();
  }

protected:
  void setSize(Real width, Real height)
  {
    mWidth  = width;
    mHeight = height;
    update();
  }

private:
  void update()
  {
    mDeltaX = mWidth / mCols;
    mDeltaY = mHeight / mRows;
    mDeltaKX = Math::twoPi / mWidth;
    mDeltaKY = Math::twoPi / mHeight;
    mKXNyquist = mDeltaKX * mCols / 2;
    mKYNyquist = mDeltaKY * mRows / 2;
  }

  int64_t mCols = 0;
  int64_t mRows = 0;
  Real    mWidth = 0;
  Real    mHeight = 0;

  // Derived quantities
  Real    mDeltaX = 0;
  Real    mDeltaY = 0;
  Real    mDeltaKX = 0;
  Real    mDeltaKY = 0;
  Real    mKXNyquist = 0;
  Real    mKYNyquist = 0;
};

class AURORA_API Buffer2DInfo : public Buffer2DInfoBase, public Energy
{
public:
  Buffer2DInfo() = default;

  Buffer2DInfo(int64_t cols, int64_t rows, Real width, Real height,
               Real energy) :
    Buffer2DInfoBase(cols, rows, width, height), Energy(energy)
  {}

  bool operator==(const Buffer2DInfo& other) const
  {
    return Buffer2DInfoBase::operator==(other) &&
           Energy::operator==(other);
  }

  bool operator!=(const Buffer2DInfo& other) const
  {
    return !operator==(other);
  }
};

} // namespace Aurora

#endif
