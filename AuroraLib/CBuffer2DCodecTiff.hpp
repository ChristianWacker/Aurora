//--- Aurora/AuroraLib/CBuffer2DCodecTiff.cpp ----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_C_BUFFER_2D_CODEC_TIFF_HPP
#define AURORA_C_BUFFER_2D_CODEC_TIFF_HPP

#include "AuroraLib/Loader.hpp"
#include "AuroraLib/Saver.hpp"

#include <string>

namespace Aurora
{

/// Class to load and save a ComplexBuffer2D in a TIFF file. The data is stored
/// in two layers with single precision (32-bit) floats. The first layer
/// represents the real part and the second layer the imaginary part of the
/// complex values.
class CBuffer2DCodecTiff :
  public Loader<CBuffer2D, Buffer2DInfo&, std::string&>,
  public Saver<CBuffer2D, const Buffer2DInfo&, const std::string&>
{
public:
  void load(const std::string& filename, CBuffer2D& buffer, Buffer2DInfo& info,
            std::string& description) override;

  void save(const CBuffer2D& buffer, const Buffer2DInfo& info,
            const std::string& description,
            const std::string& filename) override;
};

}

#endif
