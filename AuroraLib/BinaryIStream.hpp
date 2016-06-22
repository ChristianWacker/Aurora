//--- Aurora/AuroraLib/BinaryIStream.hpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_BINARY_I_STREAM_HPP
#define AURORA_AURORA_LIB_BINARY_I_STREAM_HPP

#include "AuroraLib/MemoryBlock.hpp"
#include "AuroraLib/Utils.hpp"

#include <cstring>

namespace Aurora
{

class BinaryIStream
{
public:
  static BinaryIStream openFile(const std::string& filename)
  {
    return BinaryIStream(MemoryBlock::openFile(filename));
  }

  BinaryIStream(const MemoryBlock& block) :
    mBlock(block), mPos(mBlock.data()), mEnd(mBlock.data() + mBlock.size())
  {}

  template<class T>
  T read()
  {
    T result;
    read<T>(&result, 1);
    return result;
  }

  template<class T>
  void read(T* dst, size_t count)
  {
    read(dst, sizeof(T), count);
  }

  bool checkSignature(const char* signature, size_t length = 0)
  {
    if (!length)
      length = strlen(signature);

    bool result = (0 == memcmp(mPos, signature, length));

    if (result)
      seek(length);

    return result;
  }

  void seek(intptr_t count)
  {
    mPos = clamp(mPos + count, begin(), end());
  }

  const char* begin() const { return mBlock.data(); }
  const char* end() const { return mEnd; }
  const char* pos() const { return mPos; }
  size_t size() const { return mBlock.size(); }

private:
  void read(void* buffer, size_t size, size_t count)
  {
    size_t numBytes = size * count;

    if (mPos + numBytes > mEnd)
      numBytes = mEnd - mPos;
    if (0 == numBytes)
      return;

    memcpy(buffer, mPos, numBytes);
    mPos += numBytes;
  }

  MemoryBlock mBlock;

  const char* mPos;
  const char* mEnd;
};

} // namespace Aurora

#endif
