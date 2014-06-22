//******************************************************************************
/*!
  \file      src/RNGTest/TestStack.h
  \author    J. Bakosi
  \date      Sat 14 Jun 2014 02:34:40 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
