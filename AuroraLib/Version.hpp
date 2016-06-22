//--- Aurora/AuroraLib/Version.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_VERSION_HPP
#define AURORA_AURORA_LIB_VERSION_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <string>

namespace Aurora
{

/// Aurora version.
AURORA_API extern const std::string version;

AURORA_API std::string buildString();

} // namespace Aurora

#endif
