//--- Aurora/AuroraLib/Laplace.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_LAPLACE_HPP
#define AURORA_AURORA_LIB_LAPLACE_HPP

#include "AuroraLib/Bandlimit.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/MultisliceParameters.hpp"

namespace Aurora
{

/// Types for the different point rules
struct ThreePoint { static constexpr int64_t margin = 1; };
struct FivePoint { static constexpr int64_t margin = 2; };
struct SevenPoint { static constexpr int64_t margin = 3; };
struct NinePoint { static constexpr int64_t margin = 4; };

template<class Method>
class Laplace;

template<class Method>
class LaplaceBase : public BufferBase<Laplace<Method>>
{
public:
  typedef Complex Element;

  LaplaceBase(const CMarginBuffer2D<Method::margin>& buffer) :
    mBuffer(buffer) { }

  int64_t cols() const { return mBuffer.cols(); }
  int64_t rows() const { return mBuffer.rows(); }

protected:
  const CMarginBuffer2D<Method::margin>& mBuffer;
};

template<>
class Laplace<ThreePoint> : public LaplaceBase<ThreePoint>
{
public:
  // use the base class constructor
  using LaplaceBase<ThreePoint>::LaplaceBase;

  Complex pixel(int64_t x, int64_t y) const
  {
    return   factor[0] *    mBuffer.pixel(x    , y    )
           + factor[1] * (  mBuffer.pixel(x + 1, y    )
                          + mBuffer.pixel(x - 1, y    )
                          + mBuffer.pixel(x    , y + 1)
                          + mBuffer.pixel(x    , y - 1));
  }

  #if (AURORA_COMPILER == AURORA_COMPILER_MSVC)
  static const Real factor[2];
  #else
  static constexpr Real factor[2] = {-4, 1};
  #endif
};

template<>
class Laplace<FivePoint> : public LaplaceBase<FivePoint>
{
public:
  // use the base class constructor
  using LaplaceBase<FivePoint>::LaplaceBase;

  Complex pixel(int64_t x, int64_t y) const
  {
    return   factor[0] *    mBuffer.pixel(x    , y    )
           + factor[1] * (  mBuffer.pixel(x + 1, y    )
                          + mBuffer.pixel(x - 1, y    )
                          + mBuffer.pixel(x    , y + 1)
                          + mBuffer.pixel(x    , y - 1))
           + factor[2] * (  mBuffer.pixel(x + 2, y    )
                          + mBuffer.pixel(x - 2, y    )
                          + mBuffer.pixel(x    , y + 2)
                          + mBuffer.pixel(x    , y - 2));
  }

  #if (AURORA_COMPILER == AURORA_COMPILER_MSVC)
  static const Real factor[3];
  #else
  static constexpr Real factor[3] = {-5, R(4.0) / R(3.0), R(-1.0) / R(12.0)};
  #endif
};

template<>
class Laplace<SevenPoint> : public LaplaceBase<SevenPoint>
{
public:
  // use the base class constructor
  using LaplaceBase<SevenPoint>::LaplaceBase;

  Complex pixel(int64_t x, int64_t y) const
  {
    return   factor[0] *    mBuffer.pixel(x    , y    )
           + factor[1] * (  mBuffer.pixel(x + 1, y    )
                          + mBuffer.pixel(x - 1, y    )
                          + mBuffer.pixel(x    , y + 1)
                          + mBuffer.pixel(x    , y - 1))
           + factor[2] * (  mBuffer.pixel(x + 2, y    )
                          + mBuffer.pixel(x - 2, y    )
                          + mBuffer.pixel(x    , y + 2)
                          + mBuffer.pixel(x    , y - 2))
           + factor[3] * (  mBuffer.pixel(x + 3, y    )
                          + mBuffer.pixel(x - 3, y    )
                          + mBuffer.pixel(x    , y + 3)
                          + mBuffer.pixel(x    , y - 3));
  }

  #if (AURORA_COMPILER == AURORA_COMPILER_MSVC)
  static const Real factor[4];
  #else
  static constexpr Real factor[4] =
    {R(-49.0) / R(9.0), R(1.5), R(-0.15), R(1.0) / R(90.0)};
  #endif
};

template<>
class Laplace<NinePoint> : public LaplaceBase<NinePoint>
{
public:
  // use the base class constructor
  using LaplaceBase<NinePoint>::LaplaceBase;

  Complex pixel(int64_t x, int64_t y) const
  {
    return   factor[0] *    mBuffer.pixel(x    , y    )
           + factor[1] * (  mBuffer.pixel(x + 1, y    )
                          + mBuffer.pixel(x - 1, y    )
                          + mBuffer.pixel(x    , y + 1)
                          + mBuffer.pixel(x    , y - 1))
           + factor[2] * (  mBuffer.pixel(x + 2, y    )
                          + mBuffer.pixel(x - 2, y    )
                          + mBuffer.pixel(x    , y + 2)
                          + mBuffer.pixel(x    , y - 2))
           + factor[3] * (  mBuffer.pixel(x + 3, y    )
                          + mBuffer.pixel(x - 3, y    )
                          + mBuffer.pixel(x    , y + 3)
                          + mBuffer.pixel(x    , y - 3))
           + factor[4] * (  mBuffer.pixel(x + 4, y    )
                          + mBuffer.pixel(x - 4, y    )
                          + mBuffer.pixel(x    , y + 4)
                          + mBuffer.pixel(x    , y - 4));
  }

  #if (AURORA_COMPILER == AURORA_COMPILER_MSVC)
  static const Real factor[5];
  #else
  static constexpr Real factor[5] =
    {R(-205.0) / R(36.0), R(1.6), R(-0.2), R(8.0) / R(315.0),
     R(-1.0) / R(560.0)};
  #endif
};

class LaplaceFT
{
public:
  LaplaceFT(const MultisliceParameters& params, const CBuffer2D& inBuffer,
            CBuffer2D& outBuffer, Complex prefactor) :
    mLaplaceBuffer(params.cols(), params.rows()),
    mBufferFft(params.cols(), params.rows()),
    mFft(Fft::create(inBuffer, mBufferFft, Direction::forward, false)),
    mFftInv(Fft::create(mBufferFft, outBuffer, Direction::backward, false))
  {
    const int64_t cols = params.cols();
    const int64_t rows = params.rows();

    const Real deltaKX = params.deltaKX();
    const Real deltaKY = params.deltaKY();

    for (int64_t j = 0; j < rows; ++j)
    {
      // coordinates in reciprocal space
      const Real kY = deltaKY * (j > rows / 2 ? j - rows : j);

      for (int64_t i = 0; i < cols; ++i)
      {
        const Real kX = deltaKX * (i > cols / 2 ? i - cols : i);
        // spatial frequency squared
        const Real kPerpSqr = kX * kX + kY * kY;
        mLaplaceBuffer.pixel(i, j) = prefactor * kPerpSqr;
      }
    }

    if (Bandlimit::Disabled == params.bandlimit())
    {
      // Correct the FFT scaling
      mLaplaceBuffer *= 1 / params.numPixels<Real>();
    }
    else
    {
      // implement the bandlimit. We do not need to explicitly incorporate the
      // FFT scaling factor, as it is already included in the masks.
      RBuffer2D mask;
      if (Bandlimit::Sharp == params.bandlimit())
        mask = CBandlimitSharp::mask(params);
      else if (Bandlimit::Smooth == params.bandlimit())
        mask = CBandlimitSmooth::mask(params);
      else
        AURORA_UNREACHABLE;

      mLaplaceBuffer *= mask;
    }
  }

  void apply()
  {
    mFft->transform();
    mBufferFft.assign(mLaplaceBuffer * mBufferFft);
    mFftInv->transform();
  }

private:
  CBuffer2D mLaplaceBuffer;
  CBuffer2D mBufferFft;
  const std::unique_ptr<Fft> mFft;
  const std::unique_ptr<Fft> mFftInv;
};

} // namespace Aurora

#endif
