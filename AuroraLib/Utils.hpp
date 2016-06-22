//--- Aurora/AuroraLib/Utils.hpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_UTILS_HPP
#define AURORA_AURORA_LIB_UTILS_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <algorithm>
#include <string>

namespace Aurora
{

/// Clamps a value to the range [min, max]
template<class T>
T clamp(T x, T min, T max)
{
  return std::min(std::max(x, min), max);
}

inline uint8_t saturate(int x)
{
  return std::min(std::max(x, 0), 255);
}

inline uint8_t saturate(Real x)
{
  return saturate(static_cast<int>(x));
}

/// Strips whitespaces (or other characters) from the beginning and aend of a
/// string.
/// @details
///  This function returns a string with whitespaces stripped from the
///  beginning and and of @p str. Without the second parameter, trim() will
///  strip these characters.
/// @param str
///  The string that will be trimmed.
/// @param charlist
///  Optionally, the characters to strip can also be specified using the
///  charlist parameter. This is simply a list of all characters that should
///  be stripped.
inline std::string trim(const std::string& str,
                        const std::string& charlist = " \t\n\v\f\r")
{
  size_t strBegin = str.find_first_not_of(charlist);
  if (std::string::npos == strBegin)
    return "";

  size_t strEnd = str.find_last_not_of(charlist);
  size_t strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}

AURORA_API std::string toString(size_t x, int width, char fill = '0');
AURORA_API std::string homeDirectory();

/// Returns the path of the temporary directory. This function checks
/// successively $TMPDIR, $TMP, $TEMP and "/tmp"
AURORA_API std::string tempDirectory();

inline std::string toUpper(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
  return s;
}

} // namespace Aurora

#endif
