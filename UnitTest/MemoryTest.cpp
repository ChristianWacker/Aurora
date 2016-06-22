//--- Aurora/UnitTest/MemoryTest.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Memory.hpp"

#include <gtest/gtest.h>

namespace Aurora
{

TEST(Memory, Allocate)
{
  void* data;

  data = Memory::alignedNew<uint8_t>(1, 8);
  EXPECT_EQ(0, reinterpret_cast<intptr_t>(data) % 8);
  Memory::alignedDelete(data);

  data = Memory::alignedNew<uint8_t>(1, 32);
  EXPECT_EQ(0, reinterpret_cast<intptr_t>(data) % 32);
  Memory::alignedDelete(data);

  data = Memory::alignedNew<uint8_t>(1, 128);
  EXPECT_EQ(0, reinterpret_cast<intptr_t>(data) % 128);
  Memory::alignedDelete(data);

  EXPECT_NO_THROW(Memory::alignedDelete(nullptr));
}

} // namespace Aurora
