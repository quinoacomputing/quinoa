//******************************************************************************
/*!
  \file      src/Main/TPLInfo.h
  \author    J. Bakosi
  \date      Wed Mar 19 08:44:08 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TPL info
  \details   TPL info
*/
//******************************************************************************
#ifndef TPLInfo_h
#define TPLInfo_h

#include <string>

#include <Print.h>

namespace tk {

#ifdef HAS_MKL
//! Echo MKL (Intel Math Kernel Library) version information
void echoMKL(const tk::Print& print, const std::string& title);
#endif

//! Echo Boost library version information
void echoBoost(const tk::Print& print, const std::string& title);

//! Echo OpenMP runtime version information
void echoOpenMP(const tk::Print& print, const std::string& title);

//! Echo Zoltan library version information
//void echoZoltan(const tk::Print& print, const std::string& title);

//! Echo Silo library version information
void echoSilo(const tk::Print& print, const std::string& title);

//! Echo HDF5 library version information
void echoHDF5(const tk::Print& print, const std::string& title);

//! Echo Zlib compression library version information
void echoZlib(const tk::Print& print, const std::string& title);

} // tk::

#endif // TPLInfo_h
