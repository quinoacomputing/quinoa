// *****************************************************************************
/*!
  \file      src/Control/StringParser.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     String parser base class definition
  \details   String parser base class definition. String parser base serves as
    a base class for various string parsers, e.g., command-line parsers. It does
    generic after-parser diagnostics.
*/
// *****************************************************************************

#include <vector>
#include <iostream>
#include <cstddef>
#include <string>

#include "StringParser.h"
#include "Exception.h"

using tk::StringParser;

StringParser::StringParser( int argc, char** argv ) : m_string()
// *****************************************************************************
//  Constructor
//! \param[in] argc Number of C-style character arrays in argv
//! \param[in] argv C-style character array of character arrays
//! \details   Convert C-style character array of character arrays to a single
//!   std::string substrings separated by spaces. Exception safety: basic
//!   guarantee: if an exception is thrown, the stream is in a valid state.
//! \author  J. Bakosi
// *****************************************************************************
{
  for (int i=1; i<argc; ++i)
    m_string += std::string(argv[i]) + ' ';
}

void
StringParser::diagnostics( const tk::Print& print,
                           const std::vector< std::string >& messages )
// *****************************************************************************
//  Echo errors and warnings accumulated during parsing
//! \param[in] print Pretty printer
//! \param[in] messages Vector of strings of errors and warnings
//! \author  J. Bakosi
// *****************************************************************************
{
  bool err = false;     // signaling whether there were any errors

  // Underline errors and warnings
  std::string underline( m_string.size(), ' ' );
  for (const auto& e : messages) {

    // decide if error or warning
    char underchar = ' ';
    if (e.find( "Error" ) != std::string::npos) { err = true; underchar = '^'; }
    else if (e.find( "Warning" ) != std::string::npos) underchar = '~';

    // underline error and warning
    if (underchar == '^' || underchar == '~') {
      auto sloc = e.find( "at 1," );
      if (sloc != std::string::npos) {  // if we have location info
        // skip "at 1,"
        sloc += 5;
        // find a dot starting from after "at 1,"
        const auto eloc = e.find_first_of( '.', sloc );
        // compute location of error by extracting it from error message
        const std::size_t errloc = std::stoul( e.substr( sloc, eloc-sloc ) ) - 1;
        // find beginning of erroneous argument
        sloc = m_string.rfind( ' ', errloc-1 );
        // special-handle the first argument with no space in front of it
        if (sloc == std::string::npos) sloc = 0; else ++sloc;
        // underline error and warning differently
        for (std::size_t i=sloc; i<errloc; ++i) underline[i] = underchar;
      }
    }
  }

  // Output errors and warnings underlined (if any) to quiet stream, list errors
  // and warnings, and exit
  if (!messages.empty()) {
    print % '\n';
    print % ">>> Command line parsed: '" % m_string % "'\n";
    print % ">>>                       " % underline % "\n";
    for (const auto& e : messages) print % ">>> " % e % std::endl;   // messages
    // Exit if there were any errors
    if (err) Throw( "Error(s) occurred while parsing the command line\n" );
  }
}
