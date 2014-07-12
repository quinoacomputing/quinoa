//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Boost.h
  \author    J. Bakosi
  \date      Wed Mar 19 11:07:59 2014
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Boost info
  \details   Boost info
*/
//******************************************************************************
#ifndef BoostInfo_h
#define BoostInfo_h

#include <string>

#include <Print.h>

namespace tk {

//! Echo Boost library version information
void echoBoost(const tk::Print& print, const std::string& title);

} // tk::

#endif // BoostInfo_h
