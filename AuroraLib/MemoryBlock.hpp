//--- Aurora/AuroraLib/MemoryBlock.hpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MEMORY_BLOCK_HPP
#define AURORA_AURORA_LIB_MEMORY_BLOCK_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <string>

namespace Aurora
{

/// Class to represent a block of memory.
class MemoryBlock
{
private:
  typedef std::shared_ptr<char> Data;

public:
  /// Create a memory with the content of the specified file.
  static MemoryBlock openFile(const std::string& filename);
  static MemoryBlock newMemory(size_t size)
  {
    return MemoryBlock(size, Data(new char[size], [](char* p) { delete[] p; }));
  }
  static MemoryBlock fromMemory(size_t size, char* data)
  {
    return MemoryBlock(size, Data(data, [](char*){ }));
  }

  char* data() { return mData.get(); }
  const char* data() const { return mData.get(); }

  /// Return the size of the memory block in bytes
  size_t size() const { return mSize; }

protected:
  /// Create a new MemoryBlock
  MemoryBlock(size_t size, const Data& data) :
    mSize(size), mData(data)
  {}

private:
  size_t mSize;
  Data mData;
};

} // namespace Aurora

#endif
