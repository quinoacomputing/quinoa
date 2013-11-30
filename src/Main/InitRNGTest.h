//******************************************************************************
/*!
  \file      src/Main/InitRNGTest.h
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 05:00:18 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest-specific initialization for main
  \details   RNGTest-specific initialization for main
*/
//******************************************************************************
#ifndef InitRNGTest_h
#define InitRNGTest_h

#include <Print.h>

namespace rngtest {

//! Echo Zoltan library version information
void echoZoltan(const tk::Print& print, const std::string& title);

//! Echo Silo library version information
void echoSilo(const tk::Print& print, const std::string& title);

//! Echo TPL version information
void echoTPL(const tk::Print& print);

} // rngtest::

#endif // InitRNGTest_h
