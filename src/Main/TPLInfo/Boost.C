//******************************************************************************
/*!
  \file      src/Main/TPLInfo/Boost.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:08:10 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Boost info
  \details   Boost info
*/
//******************************************************************************

#include <sstream>

#include <Config.h>
#include <TPLInfo/Boost.h>

#include <boost/version.hpp>

void tk::echoBoost( const tk::Print& print, const std::string& title )
//******************************************************************************
//  Echo Boost C++ libraries version information
//! \author  J. Bakosi
//******************************************************************************
{
  std::stringstream version;
  version << (BOOST_VERSION / 100000) << "."
          << ((BOOST_VERSION / 100) % 1000) << "."
          << (BOOST_VERSION % 100);

  print << '\n';
  print.subsection(title);
  print.item("Version", version.str());
}
