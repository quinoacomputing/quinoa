//******************************************************************************
/*!
  \file      src/Base/Printer.C
  \author    J. Bakosi
  \date      Wed Sep 11 12:50:51 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Timer
  \details   Timer
*/
//******************************************************************************

#include <string>
#include <iostream>
#include <iomanip>

#include <boost/format.hpp>

#include <Printer.h>

using namespace quinoa;

void
Printer::header(const std::string& title) const
//******************************************************************************
//  Printer base header
//! \author  J. Bakosi
//******************************************************************************
{
  using boost::format;
  using boost::io::group;

  std::cout << format("%|=80|\n") % group(std::setfill('='), "");
  std::cout << format("%|=80|\n") % title;
  std::cout << format("%|=80|\n") % group(std::setfill('='), "");
}
