//--- Aurora/AuroraLib/StringRef.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/StringRef.hpp"

#include "AuroraLib/SmallVector.hpp"

#include <algorithm>
#include <bitset>
#include <climits>
#include <cmath>
#include <cstddef>

namespace Aurora
{

/* static */ const size_t StringRef::npos;

size_t StringRef::find(char c, size_t from) const
{
  from = std::min(from, mLength);
  if (mData)
  {
    const void* p = memchr(mData + from, c, mLength - from);

    if (p)
      return static_cast<const char*>(p) - mData;
  }

  return npos;
}

size_t StringRef::findFirstOf(StringRef charList, size_t from) const
{
  std::bitset<1 << CHAR_BIT> charBits;

  for (size_t i = 0; i != charList.size(); ++i)
    charBits.set(static_cast<uchar>(charList[i]));

  for (size_t i = std::min(from, mLength); i != mLength; ++i)
    if (charBits.test(static_cast<uchar>(mData[i])))
      return i;

  return npos;
}

size_t StringRef::findFirstNotOf(StringRef charList, size_t from) const
{
  std::bitset<1 << CHAR_BIT> charBits;

  for (size_t i = 0; charList.length() != i; ++i)
    charBits.set(static_cast<uchar>(charList[i]));

  for (size_t i = std::min(from, mLength); i != mLength; ++i)
    if (!charBits.test(static_cast<uchar>(mData[i])))
      return i;

 return npos;
}

size_t StringRef::findLastNotOf(StringRef charList, size_t from) const
{
  std::bitset<1 << CHAR_BIT> charBits;

  for (size_t i = 0; charList.length() != i; ++i)
    charBits.set(static_cast<uchar>(charList[i]));

  for (size_t i = std::min(from, mLength) - 1; npos != i; --i)
    if (!charBits.test(static_cast<uchar>(mData[i])))
      return i;

  return npos;
}

StringRef StringRef::trimmed(StringRef charList) const
{
  size_t start = findFirstNotOf(charList);
  if (StringRef::npos == start)
    return StringRef();

  size_t end = findLastNotOf(charList);
  size_t count = end - start + 1;

  return subStr(start, count);
}

Real StringRef::toReal(bool* ok) const
{
  // create a null terminated buffer
  SmallVector<char, 15> buffer(mLength + 1);
  std::copy(mData, mData + mLength, buffer.get());
  buffer[mLength] = 0;

  char* end;
  double d = strtod(buffer.get(), &end);

  if (buffer.get() + mLength != end)
  {
    if (ok)
      *ok = false;

    return nan("");
  }

  if (ok)
    *ok = true;
  return static_cast<Real>(d);
}

std::string StringRef::toLower() const
{
  std::string result(mLength, char());
  for (size_t i = 0; i != mLength; ++i)
    result[i] = asciiToLower(mData[i]);

  return result;
}

std::string StringRef::toUpper() const
{
  std::string result(mLength, char());
  for (size_t i = 0; i != mLength; ++i)
    result[i] = asciiToUpper(mData[i]);

  return result;
}

} // namespace Aurora
