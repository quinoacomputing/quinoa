//******************************************************************************
/*!
  \file      src/Main/TPLInfo/OpenMP.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:13:32 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     OpenMP info
  \details   OpenMP info
*/
//******************************************************************************
#ifndef OpenMPInfo_h
#define OpenMPInfo_h

#include <string>

#include <Print.h>

namespace tk {

//! Echo OpenMP runtime version information
void echoOpenMP( const tk::Print& print, const std::string& title );

} // tk::

#endif // OpenMPInfo_h
