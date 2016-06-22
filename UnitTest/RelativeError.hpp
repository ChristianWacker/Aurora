//--- Aurora/UnitTest/RelativeError.hpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_UNIT_TEST_RELATIVE_ERROR_HPP
#define AURORA_UNIT_TEST_RELATIVE_ERROR_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <gtest/gtest.h>

namespace Aurora
{

inline testing::AssertionResult RealRelativPredFormat(const char* expr1,
                                                      const char* expr2,
                                                      const char* relErrorExpr,
                                                      Real val1,
                                                      Real val2,
                                                      Real relError)
{
  using std::fabs;

  if (R(0.0) == val2)
  {
    if (fabs(val1) <= relError)
      return testing::AssertionSuccess();

    testing::Message msg;
    msg << "The absolute error of " << expr1  << " is " << fabs(val1)
        << ", which exceeds " << relError << ", where\n"
        << expr1 << " evaluates to " << val1 << ",\n"
        << expr2 << " evaluates to " << val2 << ", and\n"
        << relErrorExpr << " evaluates to " << relError << ".";
    return testing::AssertionFailure(msg);
  }

  if (R(0.0) == val1)
  {
    if (fabs(val2) <= relError)
      return testing::AssertionSuccess();

    testing::Message msg;
    msg << "The absolute error of " << expr2  << " is " << fabs(val2)
        << ", which exceeds " << relError << ", where\n"
        << expr1 << " evaluates to " << val1 << ",\n"
        << expr2 << " evaluates to " << val2 << ", and\n"
        << relErrorExpr << " evaluates to " << relError << ".";
    return testing::AssertionFailure(msg);
  }

  const Real relError1 = fabs((val1 - val2) / val1);
  const Real relError2 = fabs((val1 - val2) / val2);

  if ((relError1 <= relError) && (relError2 <= relError))
    return testing::AssertionSuccess();

  testing::Message msg;
  msg << "The relative errors between " << expr1 << " and " << expr2
      << " are " << relError1 << " and " << relError2 << ", which exceed "
      << relError << ", where\n"
      << expr1 << " evaluates to " << val1 << ",\n"
      << expr2 << " evaluates to " << val2 << ", and\n"
      << relErrorExpr << " evaluates to " << relError << ".";
  return AssertionFailure(msg);
}

#define EXPECT_RELATIVE(val1, val2, relError)\
  EXPECT_PRED_FORMAT3(RealRelativPredFormat, val1, val2, relError)
} // namespace testing

#endif
