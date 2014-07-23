//******************************************************************************
/*!
  \file      src/Main/TPLInfo/HDF5.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:08:38 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     HDF5 info
  \details   HDF5 info
*/
//******************************************************************************

#include <sstream>

#include <H5public.h>

#include <Config.h>
#include <TPLInfo/HDF5.h>

void tk::echoHDF5( const tk::Print& print, const std::string& title )
//******************************************************************************
//  Echo HDF5 library version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  unsigned majnum, minnum, relnum;
  H5get_libversion( &majnum, &minnum, &relnum );
  version << majnum << "." << minnum << "." << relnum;

  print << '\n';
  print.subsection(title);
  print.item("Version", version.str());
}
