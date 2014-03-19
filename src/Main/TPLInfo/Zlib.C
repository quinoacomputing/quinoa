//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Zlib.C
  \author    J. Bakosi
  \date      Wed Mar 19 11:48:12 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Zlib info
  \details   Zlib info
*/
//******************************************************************************

#include <sstream>

#include <zlib.h>

#include <Config.h>
#include <TPLInfo/Zlib.h>

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
