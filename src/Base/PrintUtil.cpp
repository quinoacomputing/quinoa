// *****************************************************************************
/*!
  \file      src/Base/PrintUtil.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     String conversion utilities
  \details   Various string conversion utilities.
*/
// *****************************************************************************

#include <string>
#include <sstream>
#include <algorithm>

#include "PrintUtil.hpp"

std::string
tk::splitLines( std::string str,
                std::string indent,
                const std::string& name,
                std::size_t width )
// *****************************************************************************
//  Clean up whitespaces and format a long string into multiple lines
//! \param[in] str String to format
//! \param[in] name String to insert before string to output
//! \param[in] indent String to use as identation
//! \param[in] width Width in characters to insert newlines for output
//! \see http://stackoverflow.com/a/6892562
//! \see http://stackoverflow.com/a/8362145
//! \see https://stackoverflow.com/a/6894764
// *****************************************************************************
{
  // remove form feeds, line feeds, carriage returns, horizontal tabs,
  // vertical tabs, see http://en.cppreference.com/w/cpp/string/byte/isspace
  str.erase(
    std::remove_if( str.begin(), str.end(),
                    []( char x ){ return std::isspace( x ) && x != ' '; } ),
    str.end() );

  // remove duplicate spaces
  str.erase(
    std::unique( str.begin(), str.end(),
                 []( char a, char b ){ return a == b && a == ' '; } ),
    str.end() );

  // format str to 'width'-character-long lines with indent
  std::istringstream in(str);
  std::ostringstream os;

  os << indent << name;
  auto current = indent.size();
  std::string word;

  while (in >> word) {
    if (current + word.size() > width) {
      os << '\n' << indent;
      current = indent.size();
    }
    os << word << ' ';
    current += word.size() + 1;
  }

  return os.str();
}

std::string
tk::baselogname( const std::string& executable )
// *****************************************************************************
// Calculate base log file name
//! \param[in] executable Name of the executable
//! \return Base log file name
// *****************************************************************************
{
  return executable + "_screen.log";
}

std::string
tk::logname( const std::string& executable, int numrestart )
// *****************************************************************************
// Construct log file name
//! \param[in] executable Name of the executable
//! \param[in] numrestart Number of times restarted
//! \return The name of the log file to use
// *****************************************************************************
{
  return baselogname( executable ) +
         (numrestart ? '.' + std::to_string(numrestart) : "" );
}
