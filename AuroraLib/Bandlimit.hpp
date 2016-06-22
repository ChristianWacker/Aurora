//--- Aurora/AuroraLib/Bandlimit.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_BANDLIMIT_HPP
#define AURORA_AURORA_LIB_BANDLIMIT_HPP

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/RBuffer2D.hpp"

namespace Aurora
{

enum class Bandlimit
{
  Disabled, Smooth, Sharp
};

/// This class enforces a bandlimit on a complex buffer (@ref CBuffer2D).
class CBandlimit
{
public:
  CBandlimit(const Buffer2DInfo& info, const CBuffer2D& bufferIn,
             CBuffer2D& bufferOut) :
    mInfo(info), mBufferIn(bufferIn), mBufferOut(bufferOut)
  {}

  /// the copy-constructor is deleted, as there is no meaningful copy semantic
  CBandlimit(const CBandlimit&) = delete;
  CBandlimit& operator=(const CBandlimit&) = delete;

  virtual ~CBandlimit() {}

  virtual void apply() = 0;

  const Buffer2DInfo& info() const
  {
    return mInfo;
  }

  static std::unique_ptr<CBandlimit> create(Bandlimit bandlimit,
                                            const Buffer2DInfo& info,
                                            const CBuffer2D& bufferIn,
                                            CBuffer2D& bufferOut);

  const CBuffer2D& bufferIn() const
  {
    return mBufferIn;
  }

  CBuffer2D& bufferOut() const
  {
    return mBufferOut;
  }

private:
  Buffer2DInfo mInfo;
  const CBuffer2D& mBufferIn;
  CBuffer2D& mBufferOut;
};

class CBandlimitDisabled : public CBandlimit
{
public:
  using CBandlimit::CBandlimit;

  void apply() override
  {
    // prevent self-assignment
    if (bufferIn().data() != bufferOut().data())
      bufferOut().assign(bufferIn());
  }
};

class CBandlimitFft : public CBandlimit
{
public:
  CBandlimitFft(const Buffer2DInfo& info, const CBuffer2D& bufferIn,
                CBuffer2D& bufferOut, RBuffer2D&& mask);

  void apply() override;

private:
  const RBuffer2D mMask;
  CBuffer2D mBufferFft;
  const std::unique_ptr<Fft> mFft;
  const std::unique_ptr<Fft> mFftInv;
};

class CBandlimitSharp : public CBandlimitFft
{
public:
  CBandlimitSharp(const Buffer2DInfo& info,
                  const CBuffer2D& bufferIn, CBuffer2D& bufferOut) :
    CBandlimitFft(info, bufferIn, bufferOut, mask(info))
  {}

  static RBuffer2D mask(const Buffer2DInfo& info);
};

class CBandlimitSmooth : public CBandlimitFft
{
public:
  CBandlimitSmooth(const Buffer2DInfo& info,
                   const CBuffer2D& bufferIn, CBuffer2D& bufferOut) :
    CBandlimitFft(info, bufferIn, bufferOut, mask(info))
  {}

  static RBuffer2D mask(const Buffer2DInfo& info);
};

class RBandlimit
{
public:
  RBandlimit(Bandlimit bandlimit, const Buffer2DInfo& info);
  /// Copy-constructor creates new internal buffers to allow for multi-threading
  RBandlimit(const RBandlimit& other);

  void apply(RBuffer2D& realBuffer);

private:
  Bandlimit mBandlimit;
  CBuffer2D mTempBuffer;
  std::unique_ptr<CBandlimit> mCBandlimit;
};

} // namespace Aurora

#endif
