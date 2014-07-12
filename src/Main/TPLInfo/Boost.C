//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Boost.C
  \author    J. Bakosi
  \date      Wed Mar 19 11:35:40 2014
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Boost info
  \details   Boost info
*/
//******************************************************************************

#include <sstream>

#include <Config.h>
#include <TPLInfo/Boost.h>

#include <boost/version.hpp>

void tk::echoBoost(const tk::Print& print, const std::string& title)
//******************************************************************************
//  Echo Boost C++ libraries version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  version << (BOOST_VERSION / 100000) << "."
          << ((BOOST_VERSION / 100) % 1000) << "."
          << (BOOST_VERSION % 100);

  print.subsection(title);
  print.item("Version", version.str());
}
