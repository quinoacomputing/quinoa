//******************************************************************************
/*!
  \file      src/Main/InitQuinoa.h
  \author    J. Bakosi
  \date      Tue 17 Dec 2013 06:59:42 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa-specific initialization for main
  \details   Quinoa-specific initialization for main
*/
//******************************************************************************
#ifndef InitQuinoa_h
#define InitQuinoa_h

#include <string>

#include <Print.h>

namespace quinoa {

//! Echo Zoltan library version information
void echoZoltan(const tk::Print& print, const std::string& title);

//! Echo Silo library version information
void echoSilo(const tk::Print& print, const std::string& title);

//! Echo HDF5 library version information
void echoHDF5(const tk::Print& print, const std::string& title);

//! Echo Zlib compression library version information
void echoZlib(const tk::Print& print, const std::string& title);

//! Echo TPL version information
void echoTPL(const tk::Print& print);

} // quinoa::

#endif // InitQuinoa_h
