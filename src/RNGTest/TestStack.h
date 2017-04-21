// *****************************************************************************
/*!
  \file      src/RNGTest/TestStack.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Stack collecting all types of random number generator statistical
     tests
  \details   Stack collecting all types of statistical tests. Currently, on
    TestU01 is interfaced. More might in the future.
*/
// *****************************************************************************
#ifndef TestStack_h
#define TestStack_h

#include "TestU01Stack.h"

namespace rngtest {

//! Stack collecting all types of random RNG statistical tests
struct TestStack {
  TestStack() : TestU01() {}
  TestU01Stack TestU01;
};

} // rngtest::

#endif // TestStack_h
