//--- Aurora/AuroraLib/BinaryOStream.hpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_BINARY_O_STREAM_HPP
#define AURORA_AURORA_LIB_BINARY_O_STREAM_HPP

#include "AuroraLib/Exception.hpp"

#include <cstdio>

namespace Aurora
{

class BinaryOStream
{
public:
  BinaryOStream(const BinaryOStream&) = delete;
  BinaryOStream(BinaryOStream&&) = default;

  ~BinaryOStream()
  {
    if (mStream)
      fclose(mStream);
  }

  static BinaryOStream createFile(const std::string& filename)
  {
    return BinaryOStream(filename);
  }

  template<class T>
  void write(T x)
  {
    write(&x, 1);
  }

  template<class T>
  void write(const T* src, size_t count)
  {
    write(src, sizeof(T), count);
  }

private:
  BinaryOStream(const std::string& filename) :
    mStream(fopen(filename.c_str(), "wb"))
  {
    if (!mStream)
      AURORA_THROW(EInOut, "File could not be created.");
  }

  void write(const void* src, size_t size, size_t count)
  {
    if (fwrite(src, size, count, mStream) < count)
      AURORA_THROW(EInOut, "Write failed.");
  }

  FILE* mStream;
};

} // namespace Aurora

#endif
