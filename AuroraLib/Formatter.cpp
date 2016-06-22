//--- Aurora/AuroraLib/Formatter.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Formatter.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

#if (AURORA_PLATFORM == AURORA_PLATFORM_LINUX)
# include <unistd.h>
#elif (AURORA_PLATFORM == AURORA_PLATFORM_WINDOWS)
# include <io.h>
# include "AuroraLib/Windows.hpp"
#endif

namespace Aurora
{

std::string Byte::toString() const
{
  if (mNumBytes > 1024ll * 1024 * 1024 * 10) // greater than 10 GB
    return std::to_string(mNumBytes / (1024 * 1024 * 1024)) + " GB";
  else if (mNumBytes > 1024 * 1024 * 10) // greater than 10 MB
    return std::to_string(mNumBytes / (1024 * 1024)) + " MB";
  else if (mNumBytes > 1024 * 10) // greater than 10 kB
    return std::to_string(mNumBytes / 1024) + " kB";
  else
    return std::to_string(mNumBytes) + " B";
}

Time Time::now()
{
  return Time(std::time(nullptr));
}

std::ostream& operator<<(std::ostream& outStream, const Time& time)
{
  // convert it into a nice structure
  tm* calenderTime = std::localtime(&time.mTime);

  char oldFill  = outStream.fill('0');
  auto oldFlags = outStream.flags();
  outStream.setf(std::ios::right);
  outStream << calenderTime
            << '-' << std::setw(2) << (calenderTime->tm_mon + 1)
            << '-' << std::setw(2) << calenderTime->tm_mday
            << ' ' << std::setw(2) << calenderTime->tm_hour
            << ':' << std::setw(2) << calenderTime->tm_min
            << ':' << std::setw(2) << calenderTime->tm_sec;

  outStream.flags(oldFlags);
  outStream.fill(oldFill);

  return outStream;
}

bool isTerminal(std::ostream& s)
{
  if (&s == &std::cout)
    return isatty(1);

  if (&s == &std::cerr)
    return isatty(2);

  return false;
}

#if (AURORA_PLATFORM == AURORA_PLATFORM_LINUX)

std::ostream& operator<<(std::ostream& stream, Color color)
{
  if (isTerminal(stream))
    stream << "\033[3" << static_cast<char>(color) << 'm';
  return stream;
}

namespace Terminal
{

std::ostream& reset(std::ostream& stream)
{
  if (isTerminal(stream))
    stream << "\033[00m";
  return stream;
}

std::ostream& bold(std::ostream& stream)
{
  if (isTerminal(stream))
    stream << "\033[1m";
  return stream;
}

std::ostream& clearLine(std::ostream& stream)
{
  if (isTerminal(stream))
    stream << "\033[0G\033[2K";
  return stream;
}

} // namespace Terminal

#elif (AURORA_PLATFORM == AURORA_PLATFORM_WINDOWS)

namespace Terminal
{
  class RestoreConsole
  {
  public:
    RestoreConsole(HANDLE handle) :
      mHandle(handle), mAttributes(7)
    {
      CONSOLE_SCREEN_BUFFER_INFO info;
      if (GetConsoleScreenBufferInfo(mHandle, &info))
        mAttributes = info.wAttributes;
    }

    ~RestoreConsole()
    {
      restore();
    }

    void restore()
    {
      SetConsoleTextAttribute(mHandle, mAttributes);
    }

  private:
    HANDLE mHandle;
    WORD   mAttributes;
  };

  HANDLE getHandle(std::ostream& s)
  {
    if (&s == &std::cout)
      return GetStdHandle(STD_OUTPUT_HANDLE);

    if (&s == &std::cerr)
      return GetStdHandle(STD_ERROR_HANDLE);

    return INVALID_HANDLE_VALUE;
  }

  std::ostream& setAttribute(std::ostream& s, int mask, int value)    
  {
    if (!isTerminal(s))
      return s;

    s.flush();

    HANDLE handle = getHandle(s);
    static RestoreConsole restoreConsole(handle);

    if (-1 == value)
    {
      restoreConsole.restore();
      return s;
    }

    CONSOLE_SCREEN_BUFFER_INFO info;
    if (!GetConsoleScreenBufferInfo(handle, &info))
      return s;

    info.wAttributes &= mask;
    info.wAttributes |= value;

    SetConsoleTextAttribute(handle, info.wAttributes);
    return s;
  }

  std::ostream& reset(std::ostream& stream)
  {
    return setAttribute(stream, 0, -1);
  }

  std::ostream& bold(std::ostream& stream)
  {
    return setAttribute(stream, ~FOREGROUND_INTENSITY, FOREGROUND_INTENSITY);
  }

  std::ostream& clearLine(std::ostream& s)
  {
  #if 0
    if (!isTerminal(s))
      return s;

    s.flush();

    HANDLE handle = getHandle(s);
    CONSOLE_SCREEN_BUFFER_INFO info;
    if (!GetConsoleScreenBufferInfo(handle, &info))
      return s;

    info.dwCursorPosition.X = 0;
    SetConsoleCursorPosition(handle, info.dwCursorPosition);

    return s;
  #endif
    if (!isTerminal(s))
      return s;

    HANDLE handle = getHandle(s);
    CONSOLE_SCREEN_BUFFER_INFO info;
    if (!GetConsoleScreenBufferInfo(handle, &info))
      return s;

    s << '\r' << std::string(info.dwMaximumWindowSize.X - 1, ' ') << '\r';
    return s;
  }
}  //namespace Terminal

std::ostream& operator<<(std::ostream& stream, Color color)
{
  return Terminal::setAttribute(stream, 0xfff0, static_cast<WORD>(color));
}

#endif

} // namespace Aurora
