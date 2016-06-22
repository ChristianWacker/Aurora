//--- Aurora/AuroraLib/ImageCodecTiff.hpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_IMAGE_CODEC_TIFF_HPP
#define AURORA_AURORA_LIB_IMAGE_CODEC_TIFF_HPP

#include "AuroraLib/Loader.hpp"
#include "AuroraLib/Saver.hpp"

namespace Aurora
{

class AURORA_API ImageCodecTiff : public Loader<Image>, public Saver<Image>
{
public:
  void load(const std::string& filename, Image& image) override;
  void save(const Image& image, const std::string& filename) override;
};

}

#endif
