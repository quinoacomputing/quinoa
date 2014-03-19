//******************************************************************************
/*!
  \file      src/Main/TPLInfo/HDF5.C
  \author    J. Bakosi
  \date      Wed Mar 19 11:44:14 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     HDF5 info
  \details   HDF5 info
*/
//******************************************************************************

#include <sstream>

#include <H5public.h>

#include <Config.h>
#include <TPLInfo/HDF5.h>

void tk::echoHDF5(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo HDF5 library version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  unsigned majnum, minnum, relnum;
  H5get_libversion( &majnum, &minnum, &relnum );
  version << majnum << "." << minnum << "." << relnum;

  print.subsection(title);
  print.item("Version", version.str());
}
