//--- Aurora/AuroraLib/Windows.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_WINDOWS_HPP
#define AURORA_AURORA_LIB_WINDOWS_HPP

#include "AuroraLib/Prerequisites.hpp"

#if (AURORA_PLATFORM != AURORA_PLATFORM_WINDOWS)
#error "This file can be used on Windows only."
#endif

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#endif
