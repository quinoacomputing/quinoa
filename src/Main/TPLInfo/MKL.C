//******************************************************************************
/*!
  \file      src/Main/TPLInfo/MKL.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:08:51 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     MKL info
  \details   MKL info
*/
//******************************************************************************

#include <sstream>

#include <Config.h>
#include <TPLInfo/MKL.h>

#ifdef HAS_MKL
#include <mkl_service.h>
#endif

#ifdef HAS_MKL
void tk::echoMKL( const tk::Print& print, const std::string& title )
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

  print << '\n';
  print.subsection(title);
  print.item("Version", version.str());
  print.item("Status", vmkl.ProductStatus);
  print.item("Build", vmkl.Build);
  print.item("Platform", vmkl.Platform);
  print.item("Processor", vmkl.Processor);
}
#endif
