//--- Aurora/AuroraLib/SampleCodecAtom.hpp -------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_SAMPLE_CODEC_ATOM_HPP
#define AURORA_AURORA_LIB_SAMPLE_CODEC_ATOM_HPP

#include "AuroraLib/Loader.hpp"
#include "AuroraLib/Saver.hpp"

namespace Aurora
{

class SampleCodecAtom : public Loader<Sample>, public Saver<Sample>
{
public:
  void load(const std::string& filename, Sample& sample) override;
  void save(const Sample& sample, const std::string& filename) override;
};

} // namespace Aurora

#endif
