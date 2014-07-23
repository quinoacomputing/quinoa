//******************************************************************************
/*!
  \file      src/Main/TPLInfo/OpenMP.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:09:05 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     OpenMP info
  \details   OpenMP info
*/
//******************************************************************************

#include <sstream>

#include <Config.h>
#include <TPLInfo/OpenMP.h>

void tk::echoOpenMP( const tk::Print& print, const std::string& title )
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

  print << '\n';
  print.subsection(title);
  print.item("Version", version.str());
}
