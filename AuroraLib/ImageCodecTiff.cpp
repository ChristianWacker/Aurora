//--- Aurora/AuroraLib/ImageCodecTiff.cpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/ImageCodecTiff.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Image.hpp"
#include "AuroraLib/TiffSupport.hpp"

#include <iostream>
#include <tiffio.h>

namespace Aurora
{

void ImageCodecTiff::load(const std::string& filename, Image& image)
{
  // open the TIFF file
  Tiff::File file(TIFFOpen(filename.c_str(), "r"));

  if (!file)
    AURORA_THROW(EFileNotFound, filename);

  // check for multiple layers
  int subfileType;

  try
  {
    Tiff::getField(file, TIFFTAG_SUBFILETYPE, &subfileType);
  }
  catch (Exception&)
  {
    subfileType = 0;
  }

  bool multilayer = (subfileType & 0x2) != 0;

  if (multilayer)
    AURORA_THROW(ENotSupported, "Expected single-layer TIFF.");

  // single layer
  // extract the number of rows and columns
  uint32_t cols;
  uint32_t rows;

  Tiff::getField(file, TIFFTAG_IMAGEWIDTH , &cols);
  Tiff::getField(file, TIFFTAG_IMAGELENGTH, &rows);

  // we expect an 8-bit grey image or an float image
  uint16_t bitsPerSample;
  uint16_t samplesPerPixel;
  uint16_t sampleFormat;
  if (1 != TIFFGetField(file, TIFFTAG_SAMPLEFORMAT , &sampleFormat))
    sampleFormat = 1;

  Tiff::getField(file, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
  Tiff::getField(file, TIFFTAG_BITSPERSAMPLE  , &bitsPerSample);

  if (verbosity >= 2)
  {
    std::cout << "Image-File: " << filename << '\n';
    std::cout << "- sampleFormat: " << sampleFormat << '\n';
    std::cout << "- samplesPerPixel: " << samplesPerPixel << '\n';
    std::cout << "- bitsPerSample: " << bitsPerSample << '\n';
  }

  if ((1 == sampleFormat) && (1 == samplesPerPixel) && (8 == bitsPerSample))
  {
    // 8-bit grey image
    image = Image(cols, rows, Image::PixelFormat::Gray8);

    // read the file line by line
    for (uint32_t y = 0; y < rows; ++y)
    {
      if (-1 == TIFFReadScanline(file, image.scanline<uint8_t>(y), y))
      {
        AURORA_THROW(EInOut, "Reading of the image data from the TIFF file "
                             "failed, line: " + std::to_string(y));
      }
    }

    return;
  }

  if ((3 == sampleFormat) && (1 == samplesPerPixel) && (32 == bitsPerSample))
  {
    // float grey image
    // 8-bit grey image
    image = Image(cols, rows, Image::PixelFormat::Float);

    // read the file line by line
    for (int64_t y = 0; y < rows; ++y)
    {
      if (-1 == TIFFReadScanline(file, image.scanline<float>(y), y))
      {
        AURORA_THROW(EInOut, "Reading of the image data from the TIFF file "
                             "failed, line: " + std::to_string(y));
      }
    }

    return;
  }

  AURORA_THROW(ENotSupported, "Unsupported TIFF-Format.");
}

void ImageCodecTiff::save(const Image& image, const std::string& filename)
{
  // create a new TIFF file
  Tiff::File file(TIFFOpen(filename.c_str(), "w"));

  if (!file)
    AURORA_THROW(EInOut, "File could not be created: " + filename);

  const int64_t cols = image.cols();
  const int64_t rows = image.rows();

  uint16_t samplesPerPixel;

  // some standard tags that must be set
  Tiff::setField(file, TIFFTAG_IMAGEWIDTH   , cols);
  Tiff::setField(file, TIFFTAG_IMAGELENGTH  , rows);
  Tiff::setField(file, TIFFTAG_BITSPERSAMPLE, 8);
  Tiff::setField(file, TIFFTAG_PLANARCONFIG , PLANARCONFIG_CONTIG);
  Tiff::setField(file, TIFFTAG_ROWSPERSTRIP , rows);

  if (Image::PixelFormat::Gray8 == image.pixelFormat())
  {
    samplesPerPixel = 1;
    Tiff::setField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  }
  else if (Image::PixelFormat::RGBA8 == image.pixelFormat())
  {
    samplesPerPixel = 3;
    Tiff::setField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  }
  else if (Image::PixelFormat::Float == image.pixelFormat())
  {
    AURORA_THROW(ENotSupported, "Float images cannot be saved.");
  }
  else
  {
    AURORA_UNREACHABLE;
  }

  Tiff::setField(file, TIFFTAG_SAMPLESPERPIXEL, samplesPerPixel);

  // convert the internal format into the TIFF format
  std::vector<uint8_t> tempBuffer(cols * rows * samplesPerPixel);

  if (Image::PixelFormat::Gray8 == image.pixelFormat())
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        tempBuffer[y * cols + x] = image.pixel<Gray8>(x, y);
  }
  else if (Image::PixelFormat::RGBA8 == image.pixelFormat())
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
    {
      for (int64_t x = 0; x < cols; ++x)
      {
        tempBuffer[(y * cols + x) * 3 + 0] = image.pixel<Rgba8>(x, y).r;
        tempBuffer[(y * cols + x) * 3 + 1] = image.pixel<Rgba8>(x, y).g;
        tempBuffer[(y * cols + x) * 3 + 2] = image.pixel<Rgba8>(x, y).b;
      }
    }
  }
  else
  {
    AURORA_UNREACHABLE;
  }

  // write all the image data at once
  const int64_t numBytes = cols * rows * samplesPerPixel;
  if (-1 == TIFFWriteEncodedStrip(file, 0, tempBuffer.data(), numBytes))
    AURORA_THROW(EInOut, "Writing image data to TIFF file failed.");
}

class RegisterImageCodecTiff
{
public:
  RegisterImageCodecTiff()
  {
    auto codec = std::make_shared<ImageCodecTiff>();

    Image::registerLoader("tif", codec);
    Image::registerSaver("tif", codec);
    Image::registerLoader("tiff", codec);
    Image::registerSaver("tiff", codec);
  }
} registerImageCodecTiff;

} // namespace Aurora
