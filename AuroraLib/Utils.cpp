//--- Aurora/AuroraLib/Utils.cpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Utils.hpp"

#include "AuroraLib/Exception.hpp"

#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>

namespace Aurora {

std::string toString(size_t x, int width, char fill)
{
  std::ostringstream s;
  s.fill(fill);
  s << std::right << std::setw(width) << x;
  return s.str();
}

std::string homeDirectory()
{
  const char* homeDir = getenv("HOME");
  if (homeDir)
    return std::string(homeDir);

  const char* homeDrive = getenv("HOMEDRIVE");
  homeDir = getenv("HOMEPATH");
  if (homeDrive && homeDir)
    return std::string(homeDrive) + std::string(homeDir);

  return "";
}

std::string tempDirectory()
{
  const char* tempDir = getenv("TMPDIR");
  if (tempDir)
    return std::string(tempDir) + "/";

  tempDir = getenv("TMP");
  if (tempDir)
    return std::string(tempDir) + "/";

  tempDir = getenv("TEMP");
  if (tempDir)
    return std::string(tempDir) + "/";

  struct stat pathInfo;
  if ((0 == stat("/tmp", &pathInfo)) && (pathInfo.st_mode & S_IFDIR))
    return std::string("/tmp/");

  AURORA_THROW(EInOut, "Temp directory not found.");
}

} // namespace Aurora
