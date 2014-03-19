//******************************************************************************
/*!
  \file      src/Main/OpenMP.C
  \author    J. Bakosi
  \date      Wed Mar 19 11:14:50 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     OpenMP info
  \details   OpenMP info
*/
//******************************************************************************

#include <sstream>

#include <Config.h>
#include <TPLInfo/OpenMP.h>

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
