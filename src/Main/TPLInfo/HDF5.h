//******************************************************************************
/*!
  \file      src/Main/TPLInfo/HDF5.h
  \author    J. Bakosi
  \date      Wed Mar 19 11:43:39 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     HDF5 info
  \details   HDF5 info
*/
//******************************************************************************
#ifndef HDF5Info_h
#define HDF5Info_h

#include <string>

#include <Print.h>

namespace tk {

//! Echo HDF5 library version information
void echoHDF5(const tk::Print& print, const std::string& title);

} // tk::

#endif // HDF5Info_h
