//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Zlib.h
  \author    J. Bakosi
  \date      Wed Mar 19 11:47:27 2014
  \copyright 2005-2014, Jozsef Bakosi.
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
void echoZlib(const tk::Print& print, const std::string& title);

} // tk::

#endif // ZlibInfo_h
