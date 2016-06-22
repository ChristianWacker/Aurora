//--- Aurora/AuroraLib/MemoryTracker.hpp ---------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_MEMORY_TRACKER_HPP
#define AURORA_AURORA_LIB_MEMORY_TRACKER_HPP

#if 0

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"

#include <iostream>
#include <map>

namespace Aurora
{

class MemoryTracker
{
public:
  static void alloc(size_t size, void* p)
  {
    map()[p] = size;
    total() += size;

    if (verbosity >= 3)
    {
      std::cout << "Alloc " << Byte(size) << "\n"
                   "Total: " << Byte(total()) << '\n';
    }
  }

  static void free(void* p)
  {
    if (!p)
      return;

    auto i = map().find(p);
    if (map().end() == i)
    {
      std::cout << "Invalid free\n";
      return;
    }

    size_t size = i->second;
    total() -= size;

    if (verbosity >= 3)
    {
      std::cout << "Free " << Byte(size) << "\n"
                   "Total: " << Byte(total()) << '\n';
    }
  }

private:
  static int64_t& total()
  {
    static int64_t total_ = 0;
    return total_;
  }

  typedef std::map<void*, size_t> Map;
  static Map& map()
  {
    static Map map_;
    return map_;
  }
};

#endif

} // namespace Aurora

#endif
