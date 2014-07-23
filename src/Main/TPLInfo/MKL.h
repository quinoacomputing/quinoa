//******************************************************************************
/*!
  \file      src/Main/TPLInfo/MKL.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:13:23 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     MKL info
  \details   MKL info
*/
//******************************************************************************
#ifndef MKLInfo_h
#define MKLInfo_h

#include <string>

#include <Print.h>

namespace tk {

#ifdef HAS_MKL
//! Echo MKL (Intel Math Kernel Library) version information
void echoMKL( const tk::Print& print, const std::string& title );
#endif

} // tk::

#endif // MKLInfo_h
