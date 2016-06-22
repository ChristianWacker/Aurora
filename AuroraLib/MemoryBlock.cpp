//--- Aurora/AuroraLib/MemoryBlock.cpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/MemoryBlock.hpp"

#include "AuroraLib/Exception.hpp"
#include "AuroraLib/Resource.hpp"

#include <cstring>
#include <fcntl.h>

#if (AURORA_PLATFORM == AURORA_PLATFORM_LINUX)

#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace Aurora
{

typedef UniqueResource<struct FileDescriptorTag, int> FileDescriptor;
template<>
void FileDescriptor::release() noexcept
{
  if (-1 != mResource)
    close(mResource);
}

/* static */ MemoryBlock MemoryBlock::openFile(const std::string& filename)
{
  auto handleError = [&filename](const std::string& function)
  {
    std::string errorString = "MemoryBlock::openFile failed.\n"
                              "filename: " + filename + "\n"
                              "function: " + function + "\n"
                              "message: " + strerror(errno);

    AURORA_THROW(EInOut, errorString);
  };

  // Open the file read-only
  FileDescriptor file(open(filename.c_str(), O_RDONLY));
  if (-1 == file)
    handleError("open()");

  // Determine the file size
  const off_t filesize = lseek(file, 0, SEEK_END);
  if (off_t(-1) == filesize)
    handleError("seek()");

  // Map the file into memory
  void* data = mmap(nullptr, filesize, PROT_READ, MAP_SHARED, file, 0);
  if (MAP_FAILED == data)
    handleError("mmap()");

  // The file descriptor can be closed now. The file itself will be kept open
  // until it is unmapped.

  auto release = [filesize](char* p)
  {
    munmap(p, filesize);
  };
  return MemoryBlock(filesize, Data(static_cast<char*>(data), release));
}

} // namespace Aurora

#elif (AURORA_PLATFORM == AURORA_PLATFORM_WINDOWS)

#include "AuroraLib/Windows.hpp"

namespace Aurora
{

typedef UniqueResource<struct HandleTag, HANDLE> Handle;
template<>
void Handle::release() noexcept
{
  if (INVALID_HANDLE_VALUE != mResource)
    CloseHandle(mResource);
}

/* static */ MemoryBlock MemoryBlock::openFile(const std::string& filename)
{
  auto handleError = [&filename](const std::string& function)
  {
    // todo:: error message (GetLastError)
    std::string errorString = "MemoryBlock::openFile failed.\n"
      "filename: " + filename + "\n"
      "function: " + function;

    AURORA_THROW(EInOut, errorString);
  };

  // Open the file read-only
  Handle file(CreateFile(filename.c_str(),
              GENERIC_READ,
              FILE_SHARE_READ,
              nullptr,
              OPEN_EXISTING,
              0,
              nullptr));

  if (INVALID_HANDLE_VALUE == file)
    handleError("CreateFile()");

  LARGE_INTEGER x;
  if (0 == GetFileSizeEx(file, &x))
    handleError("GetFileSizeEx()");

  size_t filesize = x.QuadPart;

  Handle mapping(CreateFileMapping(file, nullptr, PAGE_READONLY, 0, 0, nullptr));
  if (!mapping)
    handleError("CreateFileMapping()");

  void* data = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
  if (!data)
    handleError("MapViewOfFile");

  auto release = [](char* p)
  {
    UnmapViewOfFile(p);
  };
  return MemoryBlock(filesize, Data(static_cast<char*>(data), release));
}

} // namespace Aurora

#else
# error Unknown platform.
#endif
