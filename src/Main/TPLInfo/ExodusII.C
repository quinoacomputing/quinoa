//******************************************************************************
/*!
  \file      src/Main/TPLInfo/ExodusII.C
  \author    J. Bakosi
  \date      Thu Mar 27 16:22:56 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ExodusII library API info
  \details   ExodusII library API info
*/
//******************************************************************************

#include <sstream>

#include <Config.h>
#include <TPLInfo/ExodusII.h>

#include <exodusII.h>

void tk::echoExodusII(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo ExodusII library API version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream ex_api_ver;
  ex_api_ver << EX_API_VERS;
  std::stringstream nem_file_ver;
  nem_file_ver << NEMESIS_FILE_VERSION;

  print.subsection(title);
  print.item("ExodusII API version", ex_api_ver.str());
  print.item("Nemesis file version", nem_file_ver.str());
}
