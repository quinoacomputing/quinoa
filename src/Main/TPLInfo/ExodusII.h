//******************************************************************************
/*!
  \file      src/Main/TPLInfo/ExodusII.h
  \author    J. Bakosi
  \date      Thu Mar 27 16:13:04 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ExodusII library API info
  \details   ExodusII library API info
*/
//******************************************************************************
#ifndef ExodusIIInfo_h
#define ExodusIIInfo_h

#include <string>

#include <Print.h>

namespace tk {

//! Echo ExodusII library API version information
void echoExodusII(const tk::Print& print, const std::string& title);

} // tk::

#endif // ExodusIIInfo_h
