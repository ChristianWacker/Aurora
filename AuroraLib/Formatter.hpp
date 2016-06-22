//--- Aurora/AuroraLib/Formatter.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------
//
/// @file
/// @brief This file contains functions that help by the creation of human
///  readable output.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_FORMATTER_HPP
#define AURORA_AURORA_LIB_FORMATTER_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <chrono>
#include <iomanip>
#include <ostream>
#include <sstream>

namespace Aurora
{

/// Class to convert bytes into human readable output. This class can be
/// be streamed into an ostream or converted to string with @ref toString().
class AURORA_API Byte
{
public:
  explicit Byte(size_t numBytes) :
    mNumBytes(numBytes)
  {}

  std::string toString() const;

private:
  size_t mNumBytes;
};

/// Overload of the stream operator for the class Byte.
inline std::ostream& operator<<(std::ostream& outStream, const Byte& byte)
{
  outStream << byte.toString();
  return outStream;
}

/// This class tries to generate a nice human readable output from
/// chrono::duration. It aims to get about four significant digits.
template<class Rep, class Period>
class Duration
{
public:
  explicit Duration(const std::chrono::duration<Rep, Period>& value) :
    mValue(value)
  { }

  std::string toString() const
  {
    std::ostringstream result;
    using namespace std::chrono;

    // todo: test case
    if (mValue > hours(24))
    {
      // FIXME: VS2015
      typedef std::chrono::duration<int64_t, std::ratio<24 * 3600>> days;
      auto d = duration_cast<days>(mValue);
      auto h = duration_cast<hours>(mValue - d);
      auto min = duration_cast<minutes>(mValue - d - h);

      if (d.count() > 1)
        result << d.count() << " days ";
      else
        result << "1 day ";

      result << h.count() << " h " << min.count() << " min";
    }
    else if (mValue > hours(10))
    {
      auto h = duration_cast<hours>(mValue);
      auto min = duration_cast<minutes>(mValue - h);
      result << h.count() << " h " << min.count() << " min";
    }
    else if (mValue > hours(1))
    {
      auto h = duration_cast<hours>(mValue);
      auto min = duration_cast<minutes>(mValue - h);
      auto sec = duration_cast<seconds>(mValue - h - min);
      result << h.count() << " h " << min.count() << " min "
             << sec.count() << " s";
    }
    else if (mValue > minutes(10))
    {
      auto min = duration_cast<minutes>(mValue);
      auto sec = duration_cast<seconds>(mValue - min);
      result << min.count() << " min " << sec.count() << " s";
    }
    else if (mValue > minutes(1))
    {
      typedef std::chrono::duration<int, std::ratio<1, 10>> deciseconds;
      auto min = duration_cast<minutes>(mValue);
      auto sec = duration_cast<seconds>(mValue - min);
      auto decisec = duration_cast<deciseconds>(mValue - min - sec);
      result << min.count() << " min " << sec.count()
             << "." << decisec.count() << " s";
    }
    else if (mValue > seconds(10))
    {
      typedef std::chrono::duration<int, std::ratio<1, 100>> centiseconds;
      auto sec = duration_cast<seconds>(mValue);
      auto centisec = duration_cast<centiseconds>(mValue - sec);
      result << sec.count() << "." << std::setfill('0')
             << std::setw(2) << centisec.count() << " s";
    }
    else if (mValue > seconds(1))
    {
      auto sec = duration_cast<seconds>(mValue);
      auto millisec = duration_cast<milliseconds>(mValue - sec);
      result << sec.count() << "." << std::setfill('0')
             << std::setw(3) << millisec.count() << " s";
    }
    else if (mValue > milliseconds(100))
    {
      auto millisec = duration_cast<milliseconds>(mValue);
      auto microsec = duration_cast<microseconds>(mValue - millisec);
      result << millisec.count() << "." << std::setfill('0')
             << std::setw(1) << microsec.count() / 100 << " ms";
    }
    else if (mValue > milliseconds(10))
    {
      auto millisec = duration_cast<milliseconds>(mValue);
      auto microsec = duration_cast<microseconds>(mValue - millisec);
      result << millisec.count() << "." << std::setfill('0')
             << std::setw(2) << microsec.count() / 10 << " ms";
    }
    else
    {
      result << duration_cast<microseconds>(mValue).count() << " µs";
    }

    return result.str();
  }

private:
  std::chrono::duration<Rep, Period> mValue;
};

/// Overload of the stream operator for @ref Duration.
template<class Rep, class Period>
inline std::ostream& operator<<(std::ostream& outStream,
                                const Duration<Rep, Period>& duration)
{
  outStream << duration.toString();
  return outStream;
}

/// Helper function to create an instance of @ref Duration without explicitly
/// specifying the template parameters.
template<class Rep, class Period>
Duration<Rep, Period> duration(const std::chrono::duration<Rep, Period>& d)
{
  return Duration<Rep, Period>(d);
}

class AURORA_API Time
{
public:
  explicit Time(std::time_t time) :
    mTime(time)
  { }

  static Time now();
  AURORA_API friend std::ostream& operator<<(std::ostream& outStream,
                                             const Time& time);

private:
  std::time_t mTime;
};

AURORA_API bool isTerminal(std::ostream& s);

#if (AURORA_PLATFORM == AURORA_PLATFORM_LINUX)
enum class Color : char
{
  black   = '0',
  red     = '1',
  green   = '2',
  yellow  = '3',
  blue    = '4',
  magenta = '5',
  cyan    = '6',
  white   = '7'
};

#elif (AURORA_PLATFORM == AURORA_PLATFORM_WINDOWS)
enum class Color : uint16_t
{
  black   = 0,
  blue    = 1,
  green   = 2,
  cyan    = 3,
  red     = 4,
  magenta = 5,
  yellow  = 6,
  white   = 7
};

#endif

AURORA_API std::ostream& operator<<(std::ostream& stream, Color color);

namespace Terminal
{

AURORA_API std::ostream& reset(std::ostream& stream);
AURORA_API std::ostream& bold(std::ostream& stream);
AURORA_API std::ostream& clearLine(std::ostream& stream);

} // namespace Terminal

} // namespace Aurora

#endif

