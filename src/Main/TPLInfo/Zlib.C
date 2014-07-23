//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Zlib.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:09:36 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Zlib info
  \details   Zlib info
*/
//******************************************************************************

#include <sstream>

#include <zlib.h>

#include <Config.h>
#include <TPLInfo/Zlib.h>

void tk::echoZlib( const tk::Print& print, const std::string& title )
//******************************************************************************
//  Echo Zlib compression library version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  version << zlibVersion();

  print << '\n';
  print.subsection(title);
  print.item("Version", version.str());
}
