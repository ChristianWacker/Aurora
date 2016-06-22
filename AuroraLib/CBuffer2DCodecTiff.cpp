//--- Aurora/AuroraLib/CBuffer2DCodecTiff.cpp ----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CBuffer2DCodecTiff.hpp"

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/TiffSupport.hpp"

#include <tiffio.h>
#include <vector>

namespace Aurora
{

const ttag_t customTiffTagWidth  = 65000;
const ttag_t customTiffTagHeight = 65001;
const ttag_t customTiffTagEnergy = 65002;

const uint32_t numberCustomTags = 3;
static const TIFFFieldInfo extendedTiffFieldInfo[numberCustomTags] =
{
  {customTiffTagWidth,  1, 1, TIFF_FLOAT, FIELD_CUSTOM, 0, 0, const_cast<char*>("BufferWidth") },
  {customTiffTagHeight, 1, 1, TIFF_FLOAT, FIELD_CUSTOM, 0, 0, const_cast<char*>("BufferHeight") },
  {customTiffTagEnergy, 1, 1, TIFF_FLOAT, FIELD_CUSTOM, 0, 0, const_cast<char*>("Energy") }
};

// In case we want a chain of extensions
static TIFFExtendProc parentExtender = nullptr;

static void registerCustomTiffTags(TIFF* tif)
{
  // Install the extended Tag field info
  TIFFMergeFieldInfo(tif, extendedTiffFieldInfo, numberCustomTags);

  if (parentExtender)
    (*parentExtender)(tif);
}

static void augmentLibtiffWithCustomTags()
{
  static bool firstTime = true;
  if (!firstTime)
    return;
  firstTime = false;
  parentExtender = TIFFSetTagExtender(registerCustomTiffTags);
}

void CBuffer2DCodecTiff::load(const std::string& filename,
                              CBuffer2D& buffer,
                              Buffer2DInfo& info,
                              std::string& description)
{
  augmentLibtiffWithCustomTags();

  // open the TIFF file
  Tiff::File file(TIFFOpen(filename.c_str(), "r"));

  if (!file)
    AURORA_THROW(EFileNotFound, filename);

  // check for multiple layers
  int subfileType;
  Tiff::getField(file, TIFFTAG_SUBFILETYPE, &subfileType);

  bool multilayer = (subfileType & 0x2) != 0;

  if (!multilayer)
    AURORA_THROW(ENotSupported, "Expected multi-layer TIFF.");

  // the current page
  uint16_t currentPage;
  // total number of pages
  uint16_t numTotalPages;

  Tiff::getField(file, TIFFTAG_PAGENUMBER, &currentPage, &numTotalPages);

  if ((0 != currentPage) || (2 != numTotalPages))
    AURORA_THROW(ENotSupported, "Unsupported multi-layer format.");

  // we expect an 32-bit floating-point image
  uint16_t bitsPerSample;
  uint16_t samplesPerPixel;
  uint16_t sampleFormat;

  Tiff::getField(file, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
  Tiff::getField(file, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
  Tiff::getField(file, TIFFTAG_SAMPLEFORMAT, &sampleFormat);

  if ((32 != bitsPerSample) || (1 != samplesPerPixel) || (3 != sampleFormat))
  {
    AURORA_THROW(ENotSupported,
                 "Expected one channel 32-bit floating point TIFF.");
  }

  // extract the number of columns and rows
  uint32_t cols;
  uint32_t rows;
  float width;
  float height;
  float energy;

  Tiff::getField(file, TIFFTAG_IMAGEWIDTH , &cols);
  Tiff::getField(file, TIFFTAG_IMAGELENGTH, &rows);
  Tiff::getField(file, customTiffTagWidth , &width);
  Tiff::getField(file, customTiffTagHeight, &height);
  Tiff::getField(file, customTiffTagEnergy, &energy);

  info = Buffer2DInfo(cols, rows, width, height, energy);

  char* temp;
  Tiff::getField(file, TIFFTAG_IMAGEDESCRIPTION, &temp);
  description = temp;

  // create a temporary buffers to hold the first image
  std::vector<float> realPart(cols * rows);

  // read the real data from the file line by line
  for (uint32_t y = 0; y < rows; ++y)
  {
    if (-1 == TIFFReadScanline(file, &realPart[y * cols], y))
    {
      AURORA_THROW(EInOut,
                    "Reading of real data from the TIFF file failed, line: "
                    + std::to_string(y));
    }
  }

  // go to the next directory
  if (1 != TIFFReadDirectory(file))
    AURORA_THROW(EInOut, "Reading of the second TIFF directory failed.");

  // check the second directory for consistency
  uint16_t currentPage2;
  uint16_t numTotalPages2;
  uint16_t bitsPerSample2;
  uint16_t samplesPerPixel2;
  uint16_t sampleFormat2;
  uint32_t cols2;
  uint32_t rows2;

  Tiff::getField(file, TIFFTAG_PAGENUMBER, &currentPage2,
                 &numTotalPages2);
  Tiff::getField(file, TIFFTAG_BITSPERSAMPLE, &bitsPerSample2);
  Tiff::getField(file, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel2);
  Tiff::getField(file, TIFFTAG_SAMPLEFORMAT, &sampleFormat2);
  Tiff::getField(file, TIFFTAG_IMAGEWIDTH, &cols2);
  Tiff::getField(file, TIFFTAG_IMAGELENGTH, &rows2);

  if ((2 != numTotalPages2) || (1 != currentPage2) ||
      (bitsPerSample != bitsPerSample2) ||
      (samplesPerPixel != samplesPerPixel2) ||
      (sampleFormat != sampleFormat2) || (cols != cols2) ||
      (rows != rows2))
    AURORA_THROW(EInOut, "The second directory is inconsistent.");

  buffer = CBuffer2D(cols, rows);

  // a temporary buffer to hold one line from the TIFF file
  std::vector<float> line(cols);

  // read one line of the imaginary data and store it together with the real
  // data in the pixelbuffer
  for (uint32_t y = 0; y < rows; ++y)
  {
    if (-1 == TIFFReadScanline(file, line.data(), y))
    {
      AURORA_THROW(EInOut, "Reading of imaginary data from the TIFF file "
                           "failed, line: " + std::to_string(y));
    }

    for (uint32_t x = 0; x < cols; ++x)
      buffer.pixel(x, y) = Complex(realPart[y * cols + x], line[x]);
  }
}

void CBuffer2DCodecTiff::save(const CBuffer2D& buffer,
                              const Buffer2DInfo& info,
                              const std::string& description,
                              const std::string& filename)
{
  augmentLibtiffWithCustomTags();

  if ((info.cols() != buffer.cols()) || (info.rows() != buffer.rows()))
    AURORA_THROW(EInvalidParameter, "Cols or rows do not match.");

  // todo: error checking
  // create a new object to handle a TIFF file
  Tiff::File file(TIFFOpen(filename.c_str(), "w"));
  if (!file)
    AURORA_THROW(EInOut, "File could not be created " + filename);

  const uint32_t cols = static_cast<uint32_t>(buffer.cols());
  const uint32_t rows = static_cast<uint32_t>(buffer.rows());

  // some standard tags that must be set
  Tiff::setField(file, TIFFTAG_IMAGEWIDTH, cols);
  Tiff::setField(file, TIFFTAG_IMAGELENGTH, rows);
  // grey
  Tiff::setField(file, TIFFTAG_SAMPLESPERPIXEL, 1);
  // floating-point
  Tiff::setField(file, TIFFTAG_SAMPLEFORMAT, 3);
  // single-precision
  Tiff::setField(file, TIFFTAG_BITSPERSAMPLE, 32);
  Tiff::setField(file, TIFFTAG_ROWSPERSTRIP, rows);
  Tiff::setField(file, customTiffTagWidth , static_cast<double>(info.width()));
  Tiff::setField(file, customTiffTagHeight, static_cast<double>(info.height()));
  Tiff::setField(file, customTiffTagEnergy, static_cast<double>(info.energy()));
  Tiff::setField(file, TIFFTAG_IMAGEDESCRIPTION, description.c_str());

  Tiff::setField(file, TIFFTAG_SUBFILETYPE, 2);
  Tiff::setField(file, TIFFTAG_PAGENUMBER, 0, 2);

  // we save the real and the imaginary part of the complex data separately
  // we begin with real part of the complex number
  std::vector<float> line(cols);

  for (uint32_t y = 0; y < rows; ++y)
  {
    // do the necessary format conversions
    #pragma omp parallel for
    for (int64_t x = 0; x < buffer.cols(); ++x)
      line[x] = static_cast<float>(buffer.pixel(x, y).real());

    if (-1 == TIFFWriteScanline(file, line.data(), y))
      AURORA_THROW(EInOut, "Writing real data to TIFF file failed.");
  }

  // write the first directory
  if (1 != TIFFWriteDirectory(file))
    AURORA_THROW(EInOut, "Writing directory to TIFF file failed.");

  // some standard tags that must be set
  Tiff::setField(file, TIFFTAG_IMAGEWIDTH, cols);
  Tiff::setField(file, TIFFTAG_IMAGELENGTH, rows);
  // grey
  Tiff::setField(file, TIFFTAG_SAMPLESPERPIXEL, 1);
  // floating-point
  Tiff::setField(file, TIFFTAG_SAMPLEFORMAT, 3);
  // single-precision
  Tiff::setField(file, TIFFTAG_BITSPERSAMPLE, 32);
  Tiff::setField(file, TIFFTAG_ROWSPERSTRIP, rows);

  Tiff::setField(file, TIFFTAG_SUBFILETYPE, 2);
  Tiff::setField(file, TIFFTAG_PAGENUMBER, 1, 2);

  // save the imaginary part
  for (uint32_t y = 0; y < rows; ++y)
  {
    // do the necessary format conversions
    #pragma omp parallel for
    for (int64_t x = 0; x < buffer.cols(); ++x)
      line[x] = static_cast<float>(buffer.pixel(x, y).imag());

    if (-1 == TIFFWriteScanline(file, line.data(), y))
      AURORA_THROW(EInOut, "Writing real data to TIFF file failed.");
  }

  // write the second directory
  if (1 != TIFFWriteDirectory(file))
    AURORA_THROW(EInOut, "Writing directory to TIFF file failed.");
}

class RegisterComplexBuffer2DCodecTiff
{
public:
  RegisterComplexBuffer2DCodecTiff()
  {
    auto codec = std::make_shared<CBuffer2DCodecTiff>();

    CBuffer2D::registerLoader("tif", codec);
    CBuffer2D::registerSaver("tif", codec);
    CBuffer2D::registerLoader("tiff", codec);
    CBuffer2D::registerSaver("tiff", codec);
  }
} registerComplexBuffer2DCodecTiff;

} // namespace Aurora
