//******************************************************************************
/*!
  \file      src/Base/Printer.C
  \author    J. Bakosi
  \date      Wed Sep 11 16:04:53 2013
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
  std::cout << format("%|=80|\n\n") % group(std::setfill('='), "");
}

void
Printer::section(const std::string& title) const
//******************************************************************************
//  Print section header
//! \author  J. Bakosi
//******************************************************************************
{
  using boost::format;
  using boost::io::group;

  std::cout << format(" * %1%:\n") % title;
  std::cout << format(" %1%\n") % std::string(title.size()+3,'-');
}
