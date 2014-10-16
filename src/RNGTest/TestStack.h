//******************************************************************************
/*!
  \file      src/RNGTest/TestStack.h
  \author    J. Bakosi
  \date      Sat 14 Jun 2014 02:34:40 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Stack of all types of statistical tests
  \details   Stack of all types of statistical tests
*/
//******************************************************************************
#ifndef TestStack_h
#define TestStack_h

#include <TestU01Stack.h>

namespace rngtest {

//! TestStack
struct TestStack {
  TestU01Stack TestU01;
};

} // rngtest::

#endif // TestStack_h
