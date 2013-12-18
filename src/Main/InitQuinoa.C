//******************************************************************************
/*!
  \file      src/Main/InitQuinoa.C
  \author    J. Bakosi
  \date      Tue 17 Dec 2013 06:59:34 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa-specific initialization for main
  \details   Quinoa-specific initialization for main
*/
//******************************************************************************

#include <sstream>

//#include <zoltan.h>
#include <silo.h>
#include <H5public.h>
#include <zlib.h>

#include <InitQuinoa.h>

// void quinoa::echoZoltan(const tk::Print& print, const std::string& title)
// //******************************************************************************
// //  Echo Zoltan library version information
// //! \author  J. Bakosi
// //******************************************************************************
// {
//   std::stringstream version;
//   version << ZOLTAN_VERSION_NUMBER;
// 
//   print.subsection(title);
//   print.item("Version", version.str());
// }

void quinoa::echoSilo(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo Silo library version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  version << DBVersion();

  print.subsection(title);
  print.item("Version", version.str());
}

void quinoa::echoHDF5(const tk::Print& print, const std::string& title)
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

void quinoa::echoZlib(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo Zlib compression library version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  version << zlibVersion();

  print.subsection(title);
  print.item("Version", version.str());
}

void quinoa::echoTPL(const tk::Print& print)
//******************************************************************************
//  Echo TPL version informaion for libs specific to Quinoa
//! \author  J. Bakosi
//******************************************************************************
{
  //print.raw("\n");
  //echoZoltan(print, "Zoltan library");
  print.raw("\n");
  echoSilo(print, "Silo library");
  print.raw("\n");
  echoHDF5(print, "HDF5 library");
  print.raw("\n");
  echoZlib(print, "Zlib compression library");
}
