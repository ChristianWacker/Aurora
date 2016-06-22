//--- Aurora/AuroraLib/Memory.hpp ----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MEMORY_HPP
#define AURORA_AURORA_LIB_MEMORY_HPP

#include "AuroraLib/Math.hpp"

# include <malloc.h>

namespace Aurora
{

/// Class to manage aligned memory. Aligned memory is required by certain SIMD
/// instructions. Furthermore, alignment can influence cache performance.
class AURORA_API Memory
{
public:
  /// The memory alignment that is used by @ref alloc(), if nothing else is
  /// specified.
  static const size_t defaultMemoryAlignment = 128;

  /// Allocates aligned memory.
  /// @tparam T
  ///  Type of the elements that should be allocated.
  /// @param numElements
  ///  The number of elements for which memory should be allocated.
  /// @param alignment
  ///  The requested alignment. This must be a power of two and must not greater
  ///  than 128.
  template<class T>
  static T* alignedNew(const size_t numElements,
                       const size_t alignment = defaultMemoryAlignment)
  {
    return static_cast<T*>(alloc(numElements * sizeof(T), alignment));
  }

  /// Frees the memory allocated with @ref alloc().
  /// @param ptr
  ///  The pointer to the memory block. As usual a nullptr is allowed.
  static void alignedDelete(void* ptr)
  {
    free(ptr);
  }

  static void* alloc(size_t size, size_t alignment = defaultMemoryAlignment);
  static void free(void* ptr);
};

#if (AURORA_PLATFORM == AURORA_PLATFORM_LINUX)

  inline void* Memory::alloc(size_t size, size_t alignment)
  {
    assert(alignment % sizeof(void*) == 0 &&
           "Alignment must be a multiple of the pointer size");
    assert(isPowerOfTwo(alignment / sizeof(void*)) &&
           "Alignment must be a power of two multiple of size(void*)");

    void* ptr;
    int result = posix_memalign(&ptr, alignment, size);
    if (0 != result)
      AURORA_THROW(EUnknown, "Bad alloc.");

    return ptr;
  }

  inline void Memory::free(void* ptr)
  {
    ::free(ptr);
  }

#elif (AURORA_PLATFORM == AURORA_PLATFORM_WINDOWS)

  inline void* Memory::alloc(size_t size, size_t alignment)
  {
    assert(alignment % sizeof(void*) == 0 &&
           "Alignment must be a multiple of the pointer size");
    assert(isPowerOfTwo(alignment) &&
           "Alignment must be a power of two");

    void* ptr = _aligned_malloc(size, alignment);
    if (!ptr)
      AURORA_THROW(EUnknown, "Bad alloc.");

    return ptr;
  }

  inline void Memory::free(void* ptr)
  {
    _aligned_free(ptr);
  }

#else
# error Unknown platform.
#endif

} // namespace Aurora

#endif
