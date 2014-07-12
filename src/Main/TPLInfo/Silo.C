//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Silo.C
  \author    J. Bakosi
  \date      Wed Mar 19 11:40:48 2014
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Silo info
  \details   Silo info
*/
//******************************************************************************

#include <sstream>

#include <silo.h>

#include <Config.h>
#include <TPLInfo/Silo.h>

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
