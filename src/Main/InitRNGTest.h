//******************************************************************************
/*!
  \file      src/Main/InitRNGTest.h
  \author    J. Bakosi
  \date      Wed Mar 19 08:35:25 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest-specific initialization for main
  \details   RNGTest-specific initialization for main
*/
//******************************************************************************
#ifndef InitRNGTest_h
#define InitRNGTest_h

#include <Print.h>

namespace rngtest {

//! Echo TPL version information
void echoTPL(const tk::Print& print);

} // rngtest::

#endif // InitRNGTest_h
