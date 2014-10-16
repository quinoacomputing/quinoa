//******************************************************************************
/*!
  \file      src/Main/TPLInfo/HDF5.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:13:10 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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
void echoHDF5( const tk::Print& print, const std::string& title );

} // tk::

#endif // HDF5Info_h
