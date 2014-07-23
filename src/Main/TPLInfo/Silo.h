//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Silo.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:13:40 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Silo info
  \details   Silo info
*/
//******************************************************************************
#ifndef SiloInfo_h
#define SiloInfo_h

#include <string>

#include <Print.h>

namespace tk {

//! Echo Silo library version information
void echoSilo( const tk::Print& print, const std::string& title );

} // tk::

#endif // SiloInfo_h
