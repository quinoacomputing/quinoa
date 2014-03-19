//******************************************************************************
/*!
  \file      src/Main/TPLInfo/OpenMP.h
  \author    J. Bakosi
  \date      Wed Mar 19 11:13:17 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
void echoOpenMP(const tk::Print& print, const std::string& title);

} // tk::

#endif // OpenMPInfo_h
