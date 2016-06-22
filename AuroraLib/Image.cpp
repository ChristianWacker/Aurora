//--- Aurora/AuroraLib/Image.cpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Image.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Exception.hpp"
#include "AuroraLib/RBuffer2D.hpp"

namespace Aurora
{

Image Image::fromRBuffer2D(const RBuffer2D& bufferIn, ExtractorR extract)
{
  const int64_t cols = bufferIn.cols();
  const int64_t rows = bufferIn.rows();

  Image image(cols, rows, PixelFormat::RGBA8);

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      image.pixel<Rgba8>(x, y) = extract(bufferIn.pixel(x, y));

  return image;
}

Image Image::fromCBuffer2D(const CBuffer2D& bufferIn, ExtractorCGray extract)
{
  const int64_t cols = bufferIn.cols();
  const int64_t rows = bufferIn.rows();

  Image image(cols, rows, PixelFormat::Gray8);

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      image.pixel<Gray8>(x, y) = extract(bufferIn.pixel(x, y));

  return image;
}

Image Image::fromCBuffer2D(const CBuffer2D& bufferIn, ExtractorCColor extract)
{
  const int64_t cols = bufferIn.cols();
  const int64_t rows = bufferIn.rows();

  Image image(cols, rows, PixelFormat::RGBA8);

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      image.pixel<Rgba8>(x, y) = extract(bufferIn.pixel(x, y));

  return image;
}

Image::Loaders& Image::loaders()
{
  static Loaders loaders;
  return loaders;
}

Image::Savers& Image::savers()
{
  static Savers savers;
  return savers;
}

} // namespace Aurora
