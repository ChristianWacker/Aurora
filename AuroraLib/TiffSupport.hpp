//--- Aurora/AuroraLib/TiffSupport.hpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_TIFF_SUPPORT_HPP
#define AURORA_AURORA_LIB_TIFF_SUPPORT_HPP

#include "AuroraLib/Exception.hpp"
#include "AuroraLib/Resource.hpp"

#include <cstdarg>
#include <tiffio.h>

namespace Aurora
{

namespace Tiff
{
  typedef UniqueResource<struct TiffFileTag, TIFF*> File;

  /// Wrapper function for TIFFGetField to amend the function with error
  /// checking
  inline void getField(File& tiffFile, ttag_t tag, ...)
  {
    va_list arguments;
    va_start(arguments, tag);
    int result = TIFFVGetField(tiffFile, tag, arguments);
    va_end(arguments);

    if (result != 1)
      AURORA_THROW(EInOut, "Tag " + std::to_string(tag) + " not found.");
  }

  /// Wrapper function for TIFFSetField to amend the function with error
  /// checking
  inline void setField(File& tiffFile, ttag_t tag, ...)
  {
    va_list arguments;
    va_start(arguments, tag);
    int result = TIFFVSetField(tiffFile, tag, arguments);
    va_end(arguments);

    if (result != 1)
      AURORA_THROW(EInOut, "Tag " + std::to_string(tag) + " not found.");
  }
}

template<>
inline void Tiff::File::release() noexcept
{
  if (mResource)
    TIFFClose(mResource);
}

} // namespace Aurora

#endif
