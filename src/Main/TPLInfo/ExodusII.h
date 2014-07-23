//******************************************************************************
/*!
  \file      src/Main/TPLInfo/ExodusII.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:13:00 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
void echoExodusII( const tk::Print& print, const std::string& title );

} // tk::

#endif // ExodusIIInfo_h
