//--- Aurora/UnitTest/StringRefTest.cpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/StringRef.hpp"

#include <climits>
#include <gtest/gtest.h>

namespace Aurora
{

TEST(StringRef, asciiIsDigit)
{
  EXPECT_TRUE(asciiIsDigit('0'));
  EXPECT_TRUE(asciiIsDigit('1'));
  EXPECT_TRUE(asciiIsDigit('9'));
  EXPECT_FALSE(asciiIsDigit('a'));
}

TEST(StringRef, asciiIsSpace)
{
  EXPECT_TRUE(asciiIsSpace(' '));
  EXPECT_TRUE(asciiIsSpace('\t'));
  EXPECT_TRUE(asciiIsSpace('\n'));
  EXPECT_TRUE(asciiIsSpace('\v'));
  EXPECT_TRUE(asciiIsSpace('\f'));
  EXPECT_TRUE(asciiIsSpace('\r'));
  EXPECT_FALSE(asciiIsSpace('9'));
  EXPECT_FALSE(asciiIsSpace('a'));
}

TEST(StringRef, asciiToLower)
{
  EXPECT_EQ('a', asciiToLower('a'));
  EXPECT_EQ('a', asciiToLower('A'));
  EXPECT_EQ('+', asciiToLower('+'));
}

TEST(StringRef, asciiToUpper)
{
  EXPECT_EQ('A', asciiToUpper('a'));
  EXPECT_EQ('A', asciiToUpper('A'));
  EXPECT_EQ('+', asciiToUpper('+'));
}

TEST(StringRef, Initialization)
{
  EXPECT_EQ("", StringRef());
  EXPECT_EQ("", StringRef(nullptr));
  EXPECT_EQ("", StringRef(""));
  EXPECT_EQ("Hello, World!", StringRef("Hello, World!"));
  EXPECT_EQ("Hello, World!", StringRef(std::string("Hello, World!")));
  EXPECT_EQ("Hello", StringRef("Hello, World!", 5));
}

TEST(StringRef, str)
{
  EXPECT_EQ(std::string("Hello, World!"), StringRef("Hello, World!").str());
}

TEST(StringRef, comparison)
{
  EXPECT_TRUE (StringRef("Hello").equals(StringRef("Hello")));
  EXPECT_TRUE (StringRef().equals(StringRef()));
  EXPECT_FALSE(StringRef("Hello").equals(StringRef("World")));
  EXPECT_FALSE(StringRef("Hello").equals(StringRef()));

  EXPECT_TRUE (StringRef("Hello") == StringRef("Hello"));
  EXPECT_TRUE (StringRef()        == StringRef());
  EXPECT_FALSE(StringRef("Hello") == StringRef("World"));
  EXPECT_FALSE(StringRef("Hello") == StringRef());

  EXPECT_FALSE(StringRef("Hello") != StringRef("Hello"));
  EXPECT_FALSE(StringRef()        != StringRef());
  EXPECT_TRUE (StringRef("Hello") != StringRef("World"));
  EXPECT_TRUE (StringRef("Hello") != StringRef());
}

TEST(StringRef, Length)
{
  EXPECT_EQ( 0, StringRef().length());
  EXPECT_EQ( 0, StringRef(nullptr).length());
  EXPECT_EQ( 0, StringRef("").length());
  EXPECT_EQ(13, StringRef("Hello, World!").length());

  EXPECT_EQ( 0, StringRef().size());
  EXPECT_EQ( 0, StringRef(nullptr).size());
  EXPECT_EQ( 0, StringRef("").size());
  EXPECT_EQ(13, StringRef("Hello, World!").size());

  EXPECT_TRUE (StringRef().empty());
  EXPECT_TRUE (StringRef(nullptr).empty());
  EXPECT_TRUE (StringRef("").empty());
  EXPECT_FALSE(StringRef("Hello, World!").empty());
}

TEST(StringRef, subStr)
{
  EXPECT_EQ(""       , StringRef("Hello, World!").subStr(7, 0));
  EXPECT_EQ("World"  , StringRef("Hello, World!").subStr(7, 5));
  EXPECT_EQ("World!" , StringRef("Hello, World!").subStr(7));
  EXPECT_EQ("World!" , StringRef("Hello, World!").subStr(7, 40));
  EXPECT_EQ(" World!", StringRef("Hello, World!").frontTruncated(6));
}

TEST(StringRef, Access)
{
  EXPECT_EQ('H', StringRef("Hello, World!")[0]);
  EXPECT_EQ('o', StringRef("Hello, World!")[4]);
  EXPECT_EQ('!', StringRef("Hello, World!")[12]);
}

TEST(StringRef, Find)
{
  EXPECT_EQ(5, StringRef("Hello, World!").find(','));
  EXPECT_EQ(StringRef::npos, StringRef("Hello, World!").find('a'));
  EXPECT_EQ( 2, StringRef("Hello, World!").find('l'));
  EXPECT_EQ( 3, StringRef("Hello, World!").find('l', 3));
  EXPECT_EQ(10, StringRef("Hello, World!").find('l', 4));

  EXPECT_EQ(4, StringRef("Hello, World!").findFirstOf("oW"));
  EXPECT_EQ(StringRef::npos, StringRef("Hello, World!").findFirstOf("abc"));
  EXPECT_EQ(4, StringRef("Hello, World!").findFirstOf("oW", 4));
  EXPECT_EQ(7, StringRef("Hello, World!").findFirstOf("oW", 5));

  EXPECT_EQ(3, StringRef("  Hello, World!").findFirstNotOf(" H"));

  EXPECT_EQ(4, StringRef("Hello, World!").findLastNotOf(" ,", 6));
}

TEST(StringRef, Trim)
{
  EXPECT_EQ("Hello, World!", StringRef(",.--Hello, World!--").trimmed("-.,"));
}

TEST(StringRef, Conversion)
{
  bool ok = false;
  EXPECT_EQ(0, StringRef("").toInt(&ok));
  EXPECT_FALSE(ok);
  EXPECT_EQ(0, StringRef("+").toInt(&ok));
  EXPECT_FALSE(ok);
  EXPECT_EQ(0, StringRef("-").toInt(&ok));
  EXPECT_FALSE(ok);
  EXPECT_EQ(2, StringRef("2").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(123, StringRef("0123").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(456, StringRef("456").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(-789, StringRef("-0789").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(0, StringRef("+0").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(0, StringRef("-0").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(7, StringRef("+7").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(-7, StringRef("-7").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(INT_MAX, StringRef("2147483647").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(0, StringRef("2147483648").toInt(&ok));
  EXPECT_FALSE(ok);
  EXPECT_EQ(INT_MIN, StringRef("-2147483648").toInt(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(0, StringRef("-2147483649").toInt(&ok));
  EXPECT_FALSE(ok);

  EXPECT_EQ(LLONG_MIN, StringRef("-9223372036854775808").toInt64(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(LLONG_MAX, StringRef("9223372036854775807").toInt64(&ok));
  EXPECT_TRUE(ok);
  EXPECT_EQ(ULLONG_MAX, StringRef("18446744073709551615").toUint64(&ok));
  EXPECT_TRUE(ok);

  EXPECT_EQ(2.5, StringRef("2.5").toReal());

  EXPECT_EQ("hello, world!", StringRef("Hello, World!").toLower());
  EXPECT_EQ("HELLO, WORLD!", StringRef("Hello, World!").toUpper());
}

} // namespace Aurora
