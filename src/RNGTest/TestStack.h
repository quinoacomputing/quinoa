//******************************************************************************
/*!
  \file      src/RNGTest/TestStack.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 04:29:18 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Stack collecting all types of random number generator statistical
     tests
  \details   Stack collecting all types of statistical tests. Currently, on
    TestU01 is interfaced. More might in the future.
*/
//******************************************************************************
#ifndef TestStack_h
#define TestStack_h

#include <TestU01Stack.h>

namespace rngtest {

//! Stack collecting all types of random RNG statistical tests
struct TestStack {
  TestU01Stack TestU01;
};

} // rngtest::

#endif // TestStack_h
