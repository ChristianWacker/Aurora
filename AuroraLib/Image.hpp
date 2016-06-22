//--- Aurora/AuroraLib/Image.hpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_IMAGE_HPP
#define AURORA_AURORA_LIB_IMAGE_HPP

#include "AuroraLib/Loader.hpp"
#include "AuroraLib/Math.hpp"
#include "AuroraLib/Memory.hpp"
#include "AuroraLib/Saver.hpp"
#include "AuroraLib/Utils.hpp"

#include <functional>
#include <map>

namespace Aurora
{

class Gray8
{
public:
  Gray8() { }

  explicit Gray8(uint8_t g_) :
    g(g_)
  { }

  operator uint8_t()
  {
    return g;
  }

  uint8_t g;
};

class Rgba8
{
public:
  Rgba8() { }

  explicit Rgba8(uint8_t grey) :
    r(grey), g(grey), b(grey), a(255)
  { }

  explicit Rgba8(Gray8 grey) :
    r(grey.g), g(grey.g), b(grey.g), a(255)
  { }

  Rgba8(uint8_t r_, uint8_t g_, uint8_t b_) :
    r(r_), g(g_), b(b_), a(255)
  { }

  Rgba8(uint8_t r_, uint8_t g_, uint8_t b_, uint8_t a_) :
    r(r_), g(g_), b(b_), a(a_)
  { }

  uint8_t r, g, b, a;
};

inline Vector3R hsv2Rgb(Vector3R hsv)
{
  hsv.x() = clamp<Real>(hsv.x(), R(0.0), Math::twoPi);
  int h = static_cast<int>(hsv.x() / Math::thirdPi);
  Real f = hsv.x() / Math::thirdPi - h;

  const Real p = hsv.z() * (R(1.0) - hsv.y());
  const Real q = hsv.z() * (R(1.0) - hsv.y() * f);
  const Real t = hsv.z() * (R(1.0) - hsv.y() * (R(1.0) - f));

  switch (h)
  {
  case 0:
  case 6: return Vector3R(hsv.z(), t, p);
  case 1: return Vector3R(q, hsv.z(), p);
  case 2: return Vector3R(p, hsv.z(), t);
  case 3: return Vector3R(p, q, hsv.z());
  case 4: return Vector3R(t, p, hsv.z());
  case 5: return Vector3R(hsv.z(), p, q);
  default: return Vector3R(R(0.0), R(0.0), R(0.0));
  }
}

inline Rgba8 extractReal(const Complex& z, Real minReal, Real maxReal)
{
  Real scaledValue = (z.real() - minReal) / (maxReal - minReal);
  return Rgba8(saturate(scaledValue * 255));
}

inline Rgba8 extractImag(const Complex& z, Real minImag, Real maxImag)
{
  Real scaledValue = (z.imag() - minImag) / (maxImag - minImag);
  return Rgba8(saturate(scaledValue * 255));
}

inline Rgba8 extractAbs(const Complex& z, Real minAbs, Real maxAbs)
{
  Real scaledValue = (abs(z) - minAbs) / (maxAbs - minAbs);
  return Rgba8(saturate(scaledValue * 255));
}

inline Rgba8 extractLog1pAbs(const Complex& z, Real minLog1pAbs,
                             Real maxLog1pAbs)
{
  Real scaledValue((std::log1p(abs(z)) - minLog1pAbs) /
                   (maxLog1pAbs - minLog1pAbs));
  return Rgba8(saturate(scaledValue * 255));
}

inline Rgba8 extractAbsSqr(const Complex& z, Real minAbsSqr, Real maxAbsSqr)
{
  Real scaledValue = (absSqr(z) - minAbsSqr) / (maxAbsSqr - minAbsSqr);
  return Rgba8(saturate(scaledValue * 255));
}

inline Rgba8 extractPhase(const Complex& z)
{
  // The phases are in the range (-pi, pi]. Therefore we shift an scale them
  // to [0, 255].
  Real scaledValue = R(0.5) + arg(z) / Math::twoPi;
  return Rgba8(saturate(scaledValue * 255));
}

inline Rgba8 extractPhaseColor(const Complex& z)
{
  Vector3R hsv(arg(z) + Math::pi, R(1.0), R(1.0));
  Vector3R rgb = hsv2Rgb(hsv);

  return Rgba8(saturate(rgb.x() * 255),
               saturate(rgb.y() * 255),
               saturate(rgb.z() * 255));
}

/// Class representing an image. This is an two-dimensional array of pixels
/// with different color channels. The underlying basic data type is always
/// uint8_t.
/// @remarks
///  This class is uncopyable as it copy-semantics is not clear. However, you
///  can use the function copyTo() to copy the content to another image.
class AURORA_API Image : public Loadable<Image>, public Savable<Image>
{
public:
  enum class PixelFormat
  {
    Gray8, RGBA8, Float
  };

  static unsigned bytesPerPixel(PixelFormat pixelFormat)
  {
    switch (pixelFormat)
    {
    case PixelFormat::Gray8: return 1;
    case PixelFormat::RGBA8: return 4;
    case PixelFormat::Float: return 4;
    default: AURORA_UNREACHABLE;
    }
  }

  static unsigned numChannels(PixelFormat pixelFormat)
  {
    switch (pixelFormat)
    {
    case PixelFormat::Gray8: return 1;
    case PixelFormat::RGBA8: return 4;
    case PixelFormat::Float: return 1;
    default: AURORA_UNREACHABLE;
    }
  }

  Image() :
    mCols(0), mRows(0), mPixelFormat(PixelFormat::Gray8), mData(nullptr)
  { }

  /// Creates an new image with the specified properties.
  /// @param cols
  ///  Number of pixel columns in the new image.
  /// @param rows
  ///  Number of pixel rows in the new image.
  /// @param pixelFormat
  ///  The internal @ref PixelFormat of the new image.
  Image(int64_t cols, int64_t rows, PixelFormat pixelFormat) :
    mCols(cols), mRows(rows), mPixelFormat(pixelFormat),
    mData(Memory::alignedNew<uint8_t>(numBytes()))
  { }

  Image(const std::string& filename)
  {
    load(filename);
  }

  /// The copy constructor is deleted
  Image(const Image&) = delete;

  /// Move constructor
  Image(Image&& other) noexcept :
    mCols(other.mCols), mRows(other.mRows),
    mPixelFormat(other.mPixelFormat), mData(other.mData)
  {
    // reset other
    other.mCols = 0;
    other.mRows = 0;
    other.mData = nullptr;
  }

  /// Copy assignment operator is deleted
  Image& operator=(const Image&) = delete;

  /// Move assignment operator
  Image& operator=(Image&& other) noexcept
  {
    // free our data
    if (mData)
      Memory::alignedDelete(mData);

    // copy from other
    mCols = other.mCols;
    mRows = other.mRows;
    mPixelFormat = other.mPixelFormat;
    mData = other.mData;

    // reset other
    other.mCols = 0;
    other.mRows = 0;
    other.mData = nullptr;

    return *this;
  }

  ~Image()
  {
    if (mData)
      Memory::alignedDelete(mData);
  }

  /// Returns the number of pixel columns.
  int64_t cols() const
  {
    return mCols;
  }

  /// Returns the number of pixel rows.
  int64_t rows() const
  {
    return mRows;
  }

  /// Returns the number of channels in the image.
  unsigned numChannels() const
  {
    return numChannels(mPixelFormat);
  }

  /// Returns the number of bytes needed to store the image.
  size_t numBytes() const
  {
    return mCols * mRows * bytesPerPixel(mPixelFormat);
  }

  /// Returns the pixel format of the image
  PixelFormat pixelFormat() const
  {
    return mPixelFormat;
  }

  /// Returns the pointer to the image data.
  void* dataPtr()
  {
    return mData;
  }

  /// Returns the pointer to the image data.
  const void* dataPtr() const
  {
    return mData;
  }

  template<class T>
  const T* scanline(int64_t y) const
  {
    return static_cast<T*>(mData) + y * mCols;
  }

  template<class T>
  T* scanline(int64_t y)
  {
    return static_cast<T*>(mData) + y * mCols;
  }

  int64_t index(int64_t x, int64_t y) const
  {
    return y * mCols + x;
  }

  template<class T>
  T pixel(int64_t x, int64_t y) const
  {
    return static_cast<T*>(mData)[y * mCols + x];
  }

  template<class T>
  T& pixel(int64_t x, int64_t y)
  {
    return static_cast<T*>(mData)[y * mCols + x];
  }

  typedef std::function<Rgba8 (Real)> ExtractorR;
  static Image fromRBuffer2D(const RBuffer2D& bufferIn, ExtractorR extract);

  typedef std::function<Gray8 (Complex)> ExtractorCGray;
  static Image fromCBuffer2D(const CBuffer2D& bufferIn, ExtractorCGray extract);

  typedef std::function<Rgba8 (Complex)> ExtractorCColor;
  static Image fromCBuffer2D(const CBuffer2D& bufferIn,
                             ExtractorCColor extract);

  static Loadable<Image>::Loaders& loaders();
  static Savable<Image>::Savers& savers();

private:
  int64_t mCols;
  int64_t mRows;
  PixelFormat mPixelFormat;

  void* mData;
};

} // namespace Aurora

#endif
