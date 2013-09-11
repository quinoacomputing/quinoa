//******************************************************************************
/*!
  \file      src/Base/Printer.C
  \author    J. Bakosi
  \date      Wed Sep 11 17:07:44 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Timer
  \details   Timer
*/
//******************************************************************************

#include <iomanip>

#include <Printer.h>

using namespace quinoa;

void
Printer::title(const std::string& title) const
//******************************************************************************
//  Print title
//! \author  J. Bakosi
//******************************************************************************
{
  using boost::format;
  using boost::io::group;

  std::cout << format("%|=80|\n") % group(std::setfill('='), "");
  std::cout << format("%|=80|\n") % title;
  std::cout << format("%|=80|\n") % group(std::setfill('='), "");
}
