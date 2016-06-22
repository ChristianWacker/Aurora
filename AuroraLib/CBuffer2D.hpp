//--- Aurora/AuroraLib/CBuffer2D.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_C_BUFFER_2D_HPP
#define AURORA_AURORA_LIB_C_BUFFER_2D_HPP

#include "AuroraLib/Fft.hpp"
#include "AuroraLib/Image.hpp"
#include "AuroraLib/Loader.hpp"
#include "AuroraLib/MarginBuffer2D.hpp"
#include "AuroraLib/Saver.hpp"

#include <limits>
#include <string>

namespace Aurora
{

/// Two dimensional buffer for complex values. This buffer can be loaded from
/// and saved into a file.
template<int64_t margin>
class CMarginBuffer2D :
  public MarginBuffer2D<CMarginBuffer2D<margin>, Complex, margin>,
  public Loadable<CMarginBuffer2D<margin>, Buffer2DInfo&, std::string&>,
  public Savable<CMarginBuffer2D<margin>, const Buffer2DInfo&,
                 const std::string&>
{
private:
  typedef MarginBuffer2D<CMarginBuffer2D<margin>, Complex, margin> Base;

public:
  CMarginBuffer2D() = default;

  // Use the constructors from the base class
  using MarginBuffer2D<CMarginBuffer2D<margin>, Complex, margin>::MarginBuffer2D;

  CMarginBuffer2D(const std::string& filename, Buffer2DInfo& info,
                  std::string& description)
  {
    static_assert(0 == margin,
                  "Margin must be zero to load a buffer from a file");
    this->load(filename, info, description);
  }

  CMarginBuffer2D(CMarginBuffer2D&&) noexcept = default;
  CMarginBuffer2D& operator=(CMarginBuffer2D&&) noexcept = default;

  /// delete the copy constructor
  CMarginBuffer2D(const CMarginBuffer2D&) = delete;

  /// delete the copy assignment operator
  CMarginBuffer2D& operator=(const CMarginBuffer2D&) = delete;

  /// Return the sum of the real parts of all pixels
  Real realReduce() const
  {
    const int64_t cols = this->cols();
    const int64_t rows = this->rows();
    Real result = 0;

    #pragma omp parallel for reduction (+:result)
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        result += Aurora::real(this->pixel(x, y));

    return result;
  }

  /// Return the sum of the absolute values of all pixels.
  Real absReduce() const
  {
    const int64_t cols = this->cols();
    const int64_t rows = this->rows();
    Real result = 0;

    #pragma omp parallel for reduction (+:result)
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        result += Aurora::abs(this->pixel(x, y));

    return result;
  }

  /// Return the sum of the absolute squares of all pixels.
  Real absSqrReduce() const
  {
    const int64_t cols = this->cols();
    const int64_t rows = this->rows();
    Real result = 0;

    #pragma omp parallel for reduction (+:result)
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        result += Aurora::absSqr(this->pixel(x, y));

    return result;
  }

  CMarginBuffer2D<margin>& operator*=(Real value)
  {
    const int64_t numElements = this->numElements();
    Complex* const data = this->data();

    #pragma omp parallel for
    for (int64_t i = 0; i < numElements; ++i)
      data[i] *= value;

    return *this;
  }

  using Base::operator*=;

  /// Creates a new ComplexBuffer and fills it with the absolute values of this
  /// ComplexBuffer2D
  CMarginBuffer2D<margin> bufferAbs() const
  {
    // create a buffer to hold the absolute values of the first buffer
    auto result = this->createCompatibleBuffer();

    const size_t numElements = this->numElements();
    const Complex* __restrict src = this->data();
    Complex* __restrict dst = result.data();

    // copy the absolute values
    #pragma omp parallel for
    for (size_t i = 0; i < numElements; ++i)
      dst[i] = Complex(abs(src[i]));

    return result;
  }

  /// Creates a new ComplexBuffer and fills it with the absolute square values
  /// of this ComplexBuffer2D-
  CMarginBuffer2D<margin> bufferAbsSqr() const
  {
    // create a buffer to hold the absolute values of the first buffer
    auto result = this->createCompatibleBuffer();

    const size_t numElements = this->numElements();
    const Complex* __restrict src = this->data();
    Complex* __restrict dst = result.data();

    // copy the absolute values
    #pragma omp parallel for
    for (size_t i = 0; i < numElements; ++i)
      dst[i] = Complex(Aurora::absSqr(src[i]));

    return result;
  }

  /// Returns a new buffer that contains the FFT of this buffer.
  CBuffer2D fft(Direction direction) const
  {
    static_assert(0 == margin, "Margin must be zero");
    CBuffer2D result = this->createCompatibleBuffer();
    Fft::create(*this, result, direction, true)->transform();
    return result;
  }

  /// Performs a FFT shift.
  /// @remarks
  /// After a Fourier transformation of an image the highest energies are in the
  /// corners. However, this is not the best representation for visual analysis
  /// of an image. It is better to have the highest energies in the center of
  /// the image. Therefore, the four quadrants are usually swapped along the
  /// diagonals:
  /// <pre>
  ///   A   B               D    C
  ///              -->
  ///   C   D               B    A </pre>
  CBuffer2D fftShift() const
  {
    // This code is carefully designed to take care of buffers that have an odd
    // number of pixels in horizontal or vertical direction. We assign the
    // odd dimensions to subbuffer A.
    static_assert(0 == margin, "Margin must be zero");
    auto result = this->createCompatibleBuffer();

    const int64_t cols  = this->cols();
    const int64_t rows = this->rows();

    // A -> D
    #pragma omp parallel for
    for (int64_t y = 0; y < (rows + 1) / 2; ++y)
      for (int64_t x = 0; x < (cols + 1) / 2; ++x)
        result.pixel(x + cols / 2, y + rows / 2) = this->pixel(x, y);

    // B -> C
    #pragma omp parallel for
    for (int64_t y = 0; y < (rows + 1) / 2; ++y)
      for (int64_t x = 0; x < cols / 2; ++x)
        result.pixel(x, y + rows / 2) = this->pixel(x + (cols + 1) / 2, y);

    // C -> B
    #pragma omp parallel for
    for (int64_t y = 0; y < rows / 2; ++y)
      for (int64_t x = 0; x < (cols + 1) / 2; ++x)
        result.pixel(x + cols / 2, y) = this->pixel(x, y + (rows + 1) / 2);

    // D -> A
    #pragma omp parallel for
    for (int64_t y = 0; y < rows / 2; ++y)
    {
      for (int64_t x = 0; x < cols / 2; ++x)
      {
        result.pixel(x, y) =
          this->pixel(x + (cols + 1) / 2, y + (rows + 1) / 2);
      }
    }

    return result;
  }

  struct Info
  {
    Real minReal;
    Real maxReal;
    Real minImag;
    Real maxImag;
    Real minAbsSqr;
    Real maxAbsSqr;
    Real minAbs;
    Real maxAbs;
  };

  Info info() const
  {
    const int64_t cols = this->cols();
    const int64_t rows = this->rows();

    Info info;

    // analyze image: determine the minimal and maximal values of the real part,
    // imaginary part and the squared absolute value
    info.minReal =  std::numeric_limits<Real>::infinity();
    info.maxReal = -std::numeric_limits<Real>::infinity();
    info.minImag =  std::numeric_limits<Real>::infinity();
    info.maxImag = -std::numeric_limits<Real>::infinity();
    info.minAbsSqr =  std::numeric_limits<Real>::infinity();
    info.maxAbsSqr = -std::numeric_limits<Real>::infinity();

    for (int64_t y = 0; y < rows; ++y)
    {
      for (int64_t x = 0; x < cols; ++x)
      {
        using Aurora::absSqr;
        info.minReal = std::min(info.minReal, this->pixel(x, y).real());
        info.maxReal = std::max(info.maxReal, this->pixel(x, y).real());
        info.minImag = std::min(info.minImag, this->pixel(x, y).imag());
        info.maxImag = std::max(info.maxImag, this->pixel(x, y).imag());
        info.minAbsSqr = std::min(info.minAbsSqr, absSqr(this->pixel(x, y)));
        info.maxAbsSqr = std::max(info.maxAbsSqr, absSqr(this->pixel(x, y)));
      }
    }

    info.minAbs = std::sqrt(info.minAbsSqr);
    info.maxAbs = std::sqrt(info.maxAbsSqr);

    return info;
  }

  static CMarginBuffer2D<0> fromImage(const Image& image)
  {
    static_assert(0 == margin, "Margin must be zero");
    if (Image::PixelFormat::RGBA8 == image.pixelFormat())
    {
      AURORA_THROW(ENotSupported,
                   "RGBA image cannot be loaded into a CBuffer2D");
    }

    const int64_t cols = image.cols();
    const int64_t rows = image.rows();
    CBuffer2D buffer(cols, rows);

    if (Image::PixelFormat::Float == image.pixelFormat())
    {
      #pragma omp parallel for
      for (int64_t y = 0; y < rows; ++y)
        for (int64_t x = 0; x < cols; ++x)
          buffer.pixel(x, y) = Complex(image.pixel<float>(x, y));
    }
    else if (Image::PixelFormat::Gray8 == image.pixelFormat())
    {
      #pragma omp parallel for
      for (int64_t y = 0; y < rows; ++y)
        for (int64_t x = 0; x < cols; ++x)
          buffer.pixel(x, y) = Complex(image.pixel<Gray8>(x, y) / R(255.0));
    }
    else
    {
      AURORA_UNREACHABLE;
    }

    return buffer;
  }

  typedef CMarginBuffer2D<margin> Type;
  static typename Loadable<Type, Buffer2DInfo&,
                           std::string&>::Loaders& loaders();
  static typename Savable<Type, const Buffer2DInfo&,
                          const std::string&>::Savers& savers();
};

} // namespace Aurora

#endif
