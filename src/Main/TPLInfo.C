//******************************************************************************
/*!
  \file      src/Main/TPLInfo.C
  \author    J. Bakosi
  \date      Wed Mar 19 08:45:04 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TPL info
  \details   TPL info
*/
//******************************************************************************

#include <sstream>

//#include <zoltan.h>
#include <silo.h>
#include <H5public.h>
#include <zlib.h>
#include <boost/version.hpp>

#include <Config.h>
#include <TPLInfo.h>

#ifdef HAS_MKL
#include <mkl_service.h>
#endif

#ifdef HAS_MKL
void tk::echoMKL(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo MKL (Intel Math Kernel Library) version information
//! \author  J. Bakosi
//******************************************************************************
{
  MKLVersion vmkl;
  mkl_get_version( &vmkl );

  std::stringstream version;
  version << vmkl.MajorVersion << "."
          << vmkl.MinorVersion << "."
          << vmkl.UpdateVersion;

  print.subsection(title);
  print.item("Version", version.str());
  print.item("Status", vmkl.ProductStatus);
  print.item("Build", vmkl.Build);
  print.item("Platform", vmkl.Platform);
  print.item("Processor", vmkl.Processor);
}
#endif

void tk::echoBoost(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo Boost C++ libraries version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  version << (BOOST_VERSION / 100000) << "."
          << ((BOOST_VERSION / 100) % 1000) << "."
          << (BOOST_VERSION % 100);

  print.subsection(title);
  print.item("Version", version.str());
}

void tk::echoOpenMP(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo OpenMP runtime version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  #ifdef _OPENMP
  version << _OPENMP;
  #else
  version << "n/a";
  #endif

  print.subsection(title);
  print.item("Version", version.str());
}

// void tk::echoZoltan(const tk::Print& print, const std::string& title)
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

void tk::echoSilo(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo Silo library version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  version << std::string( DBVersion() );

  print.subsection(title);
  print.item("Version", version.str());
}

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

void tk::echoZlib(const tk::Print& print, const std::string& title)
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
