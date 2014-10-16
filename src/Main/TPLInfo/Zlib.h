//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Zlib.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:13:48 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Zlib info
  \details   Zlib info
*/
//******************************************************************************
#ifndef ZlibInfo_h
#define ZlibInfo_h

#include <string>

#include <Print.h>

namespace tk {

//! Echo Zlib compression library version information
void echoZlib( const tk::Print& print, const std::string& title );

} // tk::

#endif // ZlibInfo_h
