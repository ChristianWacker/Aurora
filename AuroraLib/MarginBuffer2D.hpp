//--- Aurora/AuroraLib/MarginBuffer2D.hpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MARGIN_BUFFER_2D_HPP
#define AURORA_AURORA_LIB_MARGIN_BUFFER_2D_HPP

#include "AuroraLib/Complex.hpp"
#include "AuroraLib/Memory.hpp"

namespace Aurora
{

template<class Derived>
class BufferBase
{
public:
  Derived& derived()
  {
    return *static_cast<Derived*>(this);
  }

  const Derived& derived() const
  {
    return *static_cast<const Derived*>(this);
  }
};

/// Base class for the two dimensional buffers, @ref CMarginBuffer2D and
/// @ref RBuffer2D. It contains functions for the access and manipulation
/// of a two dimensional array of data. It supports the notion of a margin
/// where out-of-range coordinates are accepted and used for wrap access (see
/// @ref replicateMargin()).
/// @remarks
///  The Buffer is non-copyable. To copy the content of the buffer the functions
///  @ref clone() and @ref assign() can be used.
template<class Derived, class Scalar, int64_t margin>
class MarginBuffer2D : public BufferBase<Derived>
{
public:
  /// Default constructor
  MarginBuffer2D() : mCols(0), mRows(0), mData(nullptr) {}

  /// Constructor. Does NOT initialize the buffer
  MarginBuffer2D(int64_t cols, int64_t rows) :
    mCols(cols), mRows(rows), mData(Memory::alignedNew<Scalar>(numElements()))
  {}

  /// This constructor initializes the buffer with an initial value.
  MarginBuffer2D(int64_t cols, int64_t rows, const Scalar& initialValue) :
    mCols(cols), mRows(rows), mData(Memory::alignedNew<Scalar>(numElements()))
  {
    setValue(initialValue);
  }

  /// Move constructor
  MarginBuffer2D(MarginBuffer2D&& other) noexcept :
    mCols(other.mCols), mRows(other.mRows), mData(other.mData)
  {
    other.mCols = 0;
    other.mRows = 0;
    other.mData = nullptr;
  }

  /// Delete the copy constructor. If you want to copy the content of the
  /// buffer, use @ref assign.
  MarginBuffer2D(const MarginBuffer2D&) = delete;

  /// Destructor
  ~MarginBuffer2D()
  {
    Memory::alignedDelete(mData);
  }

  /// Move assignment operator
  MarginBuffer2D& operator=(MarginBuffer2D&& other) noexcept
  {
    // free our own memory
    Memory::alignedDelete(mData);

    // copy from other
    mCols = other.mCols;
    mRows = other.mRows;
    mData = other.mData;

    // reset other
    other.mCols = 0;
    other.mRows = 0;
    other.mData = nullptr;

    return *this;
  }

  /// Delete the copy assignment operator.
  MarginBuffer2D& operator=(const MarginBuffer2D&) = delete;

  /// Frees the memory hold by this class.
  void free()
  {
    if (mData)
    {
      mCols = 0;
      mRows = 0;
      Memory::alignedDelete(mData);
      mData = nullptr;
    }
  }

  /// Returns the columns of the buffer.
  int64_t cols() const { return mCols; }

  /// Returns the rows of the buffer.
  int64_t rows() const { return mRows; }

  /// Returns the number of pixel in the buffer excluding the margin.
  int64_t numPixels() const { return mCols * mRows; }

  /// Returns the number of Scalars in the buffer, i.e. all Scalars including
  /// the margin.
  size_t numElements() const
  {
    return (mCols + 2 * margin) * (mRows + 2 * margin);
  }

  /// Checks if the passed buffer is compatible to this buffer. Two buffers
  /// are compatible, if they have the same number of columns and rows.
  bool compatible(const Derived& other) const
  {
    return (mCols == other.mCols) && (mRows == other.mRows);
  }

  /// Sets all pixel of the buffer to a new value.
  void setValue(const Scalar& value)
  {
    const int64_t n = (mCols + 2 * margin) * (mRows + 2 * margin);
    #pragma omp parallel for
    for (int64_t i = 0; i < n; ++i)
      mData[i] = value;
  }

  void setZero()
  {
    setValue(Scalar(0));
  }

  /// Returns a pointer to the stored data.
  const Scalar* data() const
  {
    return mData;
  }

  /// Returns a pointer to the stored data.
  Scalar* data()
  {
    return mData;
  }

  Scalar pixel(int64_t x, int64_t y) const
  {
    return mData[(y + margin) * (mCols + 2 * margin) + (x + margin)];
  }

  Scalar& pixel(int64_t x, int64_t y)
  {
    return mData[(y + margin) * (mCols + 2 * margin) + (x + margin)];
  }

  void replicateMargin()
  {
    // FIXME
#   if (AURORA_COMPILER == AURORA_COMPILER_MSVC)
    // this ifdef is necessary as MSVC does not implement the two phase
    // template parsing correctly
    assert(0 != margin && "Margin must be greater zero");
#   else
    static_assert(0 != margin, "Margin must be greater zero");
#   endif

    // for the first rows
    for (int64_t dstRow = -margin; dstRow < 0; ++dstRow)
    {
      int64_t srcRow = dstRow + mRows;

      for (int64_t column = 0; column < mCols; ++column)
        pixel(column, dstRow) = pixel(column, srcRow);
    }

    // for the last rows
    for (int64_t dstRow = mRows; dstRow < mRows + margin; ++dstRow)
    {
      int64_t srcRow = dstRow - mRows;

      for (int64_t column = 0; column < mCols; ++column)
        pixel(column, dstRow) = pixel(column, srcRow);
    }

    // for all rows
    for (int64_t row = -margin; row < mRows + margin; ++row)
    {
      // replicate the first Scalars in each row
      for (int64_t dstColumn = -margin; dstColumn < 0; ++dstColumn)
      {
        int64_t srcColumn = dstColumn + mCols;
        pixel(dstColumn, row) = pixel(srcColumn, row);
      }

      // replicate the last Scalars in each row
      for (int64_t dstColumn = mCols; dstColumn < mCols + margin;
           ++dstColumn)
      {
        int64_t srcColumn = dstColumn - mCols;
        pixel(dstColumn, row) = pixel(srcColumn, row);
      }
    }
  }

  /// Creates a new buffer that is compatible with this buffer.
  Derived createCompatibleBuffer() const
  {
    return Derived(mCols, mRows);
  }

  /// Creates a new buffer that contains a copy of the data.
  Derived clone() const
  {
    auto result = createCompatibleBuffer();
    result.assign(*this);
    return result;
  }

  Derived& operator*=(const Scalar& value)
  {
    const int64_t n = (mCols + 2 * margin) * (mRows + 2 * margin);

    #pragma omp parallel for
    for (int64_t i = 0; i < n; ++i)
      mData[i] *= value;

    return this->derived();
  }

  template<class OtherDerived>
  Derived& assign(const BufferBase<OtherDerived>& rhs)
  {
    const OtherDerived& other = rhs.derived();
    assert((other.cols() == mCols) && (other.rows() == mRows) &&
           "Nodes must be compatible");

    #pragma omp parallel for
    for (int64_t y = 0; y < mRows; ++y)
      for (int64_t x = 0; x < mCols; ++x)
        pixel(x, y) = other.pixel(x, y);

    return this->derived();
  }

  template<class OtherDerived>
  Derived& operator+=(const BufferBase<OtherDerived>& rhs)
  {
    const OtherDerived& other = rhs.derived();
    assert((other.cols() == mCols) && (other.rows() == mRows) &&
           "Nodes must be compatible");

    #pragma omp parallel for
    for (int64_t y = 0; y < mRows; ++y)
      for (int64_t x = 0; x < mCols; ++x)
        pixel(x, y) += other.pixel(x, y);

    return this->derived();
  }

  template<class OtherDerived>
  Derived& operator*=(const BufferBase<OtherDerived>& rhs)
  {
    const OtherDerived& other = rhs.derived();
    assert((other.cols() == mCols) && (other.rows() == mRows) &&
           "Nodes must be compatible");

    #pragma omp parallel for
    for (int64_t y = 0; y < mRows; ++y)
      for (int64_t x = 0; x < mCols; ++x)
        pixel(x, y) *= other.pixel(x, y);

    return this->derived();
  }

private:
  int64_t mCols;
  int64_t mRows;

  Scalar* mData;
};

template<class Derived, class Other>
class BufferUnaryOp : public BufferBase<Derived>
{
public:
  BufferUnaryOp(const Other& other) :
    mOther(other)
  { }

  int64_t cols() const { return mOther.cols(); }
  int64_t rows() const { return mOther.rows(); }

protected:
  const Other& mOther;
};

template<class Derived, class Left, class Right>
class BufferBinaryOp : public BufferBase<Derived>
{
public:
  BufferBinaryOp(const Left& lhs, const Right& rhs) :
    mLhs(lhs), mRhs(rhs)
  {
    assert((lhs.cols() == rhs.cols()) && (lhs.rows() == mRhs.rows()) &&
           "Nodes must be compatible");
  }

  // Returns the size of the buffer on the right hand side of the operator.
  // This arbitrary, as the two operands must be compatible anyway.
  int64_t cols() const { return mRhs.cols(); }
  int64_t rows() const { return mRhs.rows(); }

protected:
  const Left& mLhs;
  const Right& mRhs;
};

template<class Other>
class BufferAbs : public BufferUnaryOp<BufferAbs<Other>, Other>
{
public:
  using BufferUnaryOp<BufferAbs<Other>, Other>::BufferUnaryOp;

  Real pixel(int64_t x, int64_t y) const
  {
    return abs(this->mOther.pixel(x, y));
  }
};

template<class Other>
BufferAbs<Other> bufferAbs(const BufferBase<Other>& other)
{
  return BufferAbs<Other>(other.derived());
}

template<class Other>
class BufferReal : public BufferUnaryOp<BufferReal<Other>, Other>
{
public:
  using BufferUnaryOp<BufferReal<Other>, Other>::BufferUnaryOp;

  Real pixel(int64_t x, int64_t y) const
  {
    return real(this->mOther.pixel(x, y));
  }
};

template<class Other>
BufferReal<Other> bufferReal(const BufferBase<Other>& other)
{
  return BufferReal<Other>(other.derived());
}

template<class Left, class Right>
class BufferAdd : public BufferBinaryOp<BufferAdd<Left, Right>, Left, Right>
{
public:
  using BufferBinaryOp<BufferAdd<Left, Right>, Left, Right>::BufferBinaryOp;

  auto pixel(int64_t x, int64_t y) const
    -> decltype(this->mLhs.pixel(x, y) + this->mRhs.pixel(x, y))
  {
    return this->mLhs.pixel(x, y) + this->mRhs.pixel(x, y);
  }
};

template<class Left, class Right>
BufferAdd<Left, Right> operator+(const BufferBase<Left>& lhs,
                                 const BufferBase<Right>& rhs)
{
  return BufferAdd<Left, Right>(lhs.derived(), rhs.derived());
}

template<class Left, class Right>
class BufferSub : public BufferBinaryOp<BufferSub<Left, Right>, Left, Right>
{
public:
  using BufferBinaryOp<BufferSub<Left, Right>, Left, Right>::BufferBinaryOp;

  auto pixel(int64_t x, int64_t y) const
    -> decltype(this->mLhs.pixel(x, y) - this->mRhs.pixel(x, y))
  {
    return this->mLhs.pixel(x, y) - this->mRhs.pixel(x, y);
  }
};

template<class Left, class Right>
BufferSub<Left, Right> operator-(const BufferBase<Left>& lhs,
                                 const BufferBase<Right>& rhs)
{
  return BufferSub<Left, Right>(lhs.derived(), rhs.derived());
}

template<class Left, class Right>
class BufferMul : public BufferBinaryOp<BufferMul<Left, Right>, Left, Right>
{
public:
  using BufferBinaryOp<BufferMul<Left, Right>, Left, Right>::BufferBinaryOp;

  auto pixel(int64_t x, int64_t y) const
    -> decltype(this->mLhs.pixel(x, y) * this->mRhs.pixel(x, y))
  {
    return this->mLhs.pixel(x, y) * this->mRhs.pixel(x, y);
  }
};

template<class Left, class Right>
BufferMul<Left, Right> operator*(const BufferBase<Left>& lhs,
                                 const BufferBase<Right>& rhs)
{
  return BufferMul<Left, Right>(lhs.derived(), rhs.derived());
}

template<class LeftScalar, class Right>
class BufferAddConst : public BufferBase<BufferAddConst<LeftScalar, Right>>
{
private:
  const LeftScalar& mC;
  const Right& mRhs;

public:
  BufferAddConst(const LeftScalar& c, const Right& rhs) :
    mC(c), mRhs(rhs)
  { }

  int64_t cols() const { return mRhs.cols(); }
  int64_t rows() const { return mRhs.rows(); }

  auto pixel(int64_t x, int64_t y) const
    -> decltype(mC * mRhs.pixel(x, y))
  {
    return mC + mRhs.pixel(x, y);
  }
};

template<class Right>
BufferAddConst<Real, Right> operator+(const Real& c,
                                      const BufferBase<Right>& rhs)
{
  return BufferAddConst<Real, Right>(c, rhs.derived());
}

template<class LeftScalar, class Right>
class BufferSubConst : public BufferBase<BufferSubConst<LeftScalar, Right>>
{
private:
  const LeftScalar& mC;
  const Right& mRhs;

public:
  BufferSubConst(const LeftScalar& c, const Right& rhs) :
    mC(c), mRhs(rhs)
  { }

  int64_t cols() const { return mRhs.cols(); }
  int64_t rows() const { return mRhs.rows(); }

  auto pixel(int64_t x, int64_t y) const
    -> decltype(mC * mRhs.pixel(x, y))
  {
    return mC - mRhs.pixel(x, y);
  }
};

template<class Right>
BufferSubConst<Real, Right> operator-(const Real& c,
                                      const BufferBase<Right>& rhs)
{
  return BufferSubConst<Real, Right>(c, rhs.derived());
}

/// This class is part of the expression template infrastructure. It represents
/// the multiplication of an object with a numerical constant.
template<class LeftScalar, class Right>
class BufferMulConst : public BufferBase<BufferMulConst<LeftScalar, Right>>
{
private:
  const LeftScalar& mC;
  const Right& mRhs;

public:
  BufferMulConst(const LeftScalar& c, const Right& rhs) :
    mC(c), mRhs(rhs)
  { }

  int64_t cols() const { return mRhs.cols(); }
  int64_t rows() const { return mRhs.rows(); }

  auto pixel(int64_t x, int64_t y) const
    -> decltype(mC * mRhs.pixel(x, y))
  {
    return mC * mRhs.pixel(x, y);
  }
};

template<class Right>
BufferMulConst<Real, Right> operator*(const Real& c,
                                      const BufferBase<Right>& rhs)
{
  return BufferMulConst<Real, Right>(c, rhs.derived());
}

template<class Right>
BufferMulConst<Complex, Right> operator*(const Complex& c,
                                         const BufferBase<Right>& rhs)
{
  return BufferMulConst<Complex, Right>(c, rhs.derived());
}

} // namespace Aurora

#endif
