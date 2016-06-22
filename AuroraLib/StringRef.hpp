//--- Aurora/AuroraLib/StringRef.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_STRING_REF_HPP
#define AURORA_AURORA_LIB_STRING_REF_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <string>

namespace Aurora
{

inline int compareMemory(const char* a, const char* b, size_t length)
{
  if (0 == length)
    return 0;
  return memcmp(a, b, length);
}

inline bool asciiIsDigit(char c)
{
  return (c >= '0' && c <= '9');
}

inline bool asciiIsSpace(char c)
{
  return (' ' == c) || (c >= 0x9 && c <= 0xd);
}

inline char asciiToLower(char c)
{
  if (c >= 'A' && c <= 'Z')
    return c + 'a' - 'A'; // set the sixth bit

  return c;
}

inline char asciiToUpper(char c)
{
  if (c >= 'a' && c <= 'z')
    return c + 'A' - 'a';

  return c;
}

/// @brief
/// The StringRef save the reference to a string. The string can be c-string or
/// a std-string.
class AURORA_API StringRef
{
public:
  typedef size_t size_type;
  static const size_t npos = ~static_cast<size_t>(0);

  /// Create an empty StringRef.
  StringRef() :
    mData(nullptr), mLength(0)
  { }

  /// Create a StringRef from a c-string.
  /* implicit */ StringRef(const char* s) :
    mData(s), mLength(s ? strlen(s) : 0)
  { }

  StringRef(const char* s, size_t length) :
    mData(s), mLength(length)
  {
    assert((s || 0 == length) &&
           "StringRef cannot be initialized from a NULL argument with "
           "non-zero length.");
  }

  /// Create a StringRef from a std::string
  /* implicit */ StringRef(const std::string& s) :
    mData(s.data()), mLength(s.length())
  { }

  /// Convert the StringRef to a std::string.
  std::string str() const
  {
    return std::string(mData, mLength);
  }

  /// Return a pointer to the beginning of the data block containing the string.
  const char* data() const
  {
    return mData;
  }

  /// Check the strings this and @p s for equality
  bool equals(StringRef s) const
  {
    return ((mLength == s.mLength) &&
            (0 == compareMemory(mData, s.mData, mLength)));
  }

  /// Return the length of the represented string.
  size_t length() const
  {
    return mLength;
  }

  /// Return the length of the represented string.
  size_t size() const
  {
    return mLength;
  }

  /// Checks, if the length of the string is zero.
  bool empty() const
  {
    return (0 == mLength);
  }

  /// Create a StringRef for a string with the first @þ n elements dropped
  StringRef frontTruncated(size_t n = 1) const
  {
    assert(mLength >= n && "Dropping too much elements.");
    return subStr(n);
  }

  /// Create a StringRef for the substring starting at @p start with length
  /// @p count
  StringRef subStr(size_t start, size_t count = npos) const
  {
    start = std::min(start, mLength);
    return StringRef(mData + start, std::min(count, mLength - start));
  }

  /// Element-wise access
  char operator[](size_t index) const
  {
    assert(index < mLength && "Index out of range.");
    return mData[index];
  }

  /// Find the first occurence of @p c
  size_t find(char c, size_t from = 0) const;
  size_t findFirstOf(StringRef charList, size_t from = 0) const;
  size_t findFirstNotOf(StringRef charList, size_t from = 0) const;
  size_t findLastNotOf(StringRef charList, size_t from = npos) const;

  /// Create a new StringRef with @p charlist removed on the left and on the
  /// right.
  StringRef trimmed(StringRef charList = " \t\n\v\f\r") const;

  /// Convert the StringRef to a double value.
  Real toReal(bool* ok = nullptr) const;

  /// Convert the StringRef to a lower-case string.
  std::string toLower() const;

  /// Convert the StringRef to a upper-case string.
  std::string toUpper() const;

  int toInt(bool* ok = nullptr) const
  {
    int64_t n = toInt64(ok);

    if (n < std::numeric_limits<int>::min() ||
        n > std::numeric_limits<int>::max())
    {
      if (ok)
        *ok = false;
      n = 0;
    }

    return n;
  }

  int64_t toInt64(bool* ok = nullptr) const
  {
    if (empty())
    {
      if (ok)
        *ok = false;
      return 0;
    }

    uint64_t n;
    if ('+' == mData[0] || '-' == mData[0])
      n = StringRef(mData + 1, mLength - 1).toUint64(ok);
    else
      n = toUint64(ok);

    if ('-' == mData[0])
    {
      if (n > static_cast<uint64_t>(std::numeric_limits<int64_t>::max()) + 1)
      {
        if (ok)
          *ok = false;
        return std::numeric_limits<int64_t>::min();
      }
      else
      {
        return -n;
      }
    }

    if (n > std::numeric_limits<int64_t>::max())
    {
      if (ok)
        *ok = false;
      return std::numeric_limits<int64_t>::max();
    }
    else
    {
      return n;
    }
  }

  uint64_t toUint64(bool* ok = nullptr) const
  {
    if (empty())
    {
      if (ok)
        *ok = false;
      return 0;
    }

    uint64_t result = 0;
    const char* cur = mData;
    const char* end = mData + mLength;
    while (cur != end)
    {
      if (!asciiIsDigit(*cur))
      {
        if (ok)
          *ok = false;
        return 0;
      }

      result *= 10;
      result += static_cast<uint64_t>(*cur - '0');

      ++cur;
    }

    if (ok)
      *ok = true;
    return result;
  }

private:
  const char* mData;
  size_t       mLength;
};

inline bool operator==(StringRef lhs, StringRef rhs)
{
  return lhs.equals(rhs);
}

inline bool operator!=(StringRef lhs, StringRef rhs)
{
  return !lhs.equals(rhs);
}

} // namespace Aurora

#endif
