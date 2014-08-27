//******************************************************************************
/*!
  \file      src/Control/StringParser.C
  \author    J. Bakosi
  \date      Wed 27 Aug 2014 10:37:07 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     String parser
  \details   String parser
*/
//******************************************************************************

#include <vector>
#include <iostream>

#include <StringParser.h>
#include <Exception.h>

using tk::StringParser;

StringParser::StringParser(int argc, char** argv)
//******************************************************************************
//  Constructor
//! \param[in]     argc          Number of C-style character arrays in argv
//! \param[in]     argv          C-style character array of character arrays
//! \details   Convert C-style character array of character arrays to a single
//             std::string substrings separated by spaces. Exception safety:
//             basic guarantee: if an exception is thrown, the stream is in a
//             valid state.
//! \author  J. Bakosi
//******************************************************************************
{
  for (int i=1; i<argc; ++i) m_string += std::string( argv[i] ) + ' ';
}

void
StringParser::echoErrors( const tk::Print& print,
                          const std::vector< std::string >& errors )
//******************************************************************************
//  Echo errors accumulated during parsing
//! \param[in] errors Vector of strings of errors
//! \author  J. Bakosi
//******************************************************************************
{
  bool err = false;

  // Underline errors
  std::string pos( m_string.size(), ' ' );
  for (const auto& e : errors) {

    // decide if error or warning
    char underline = ' ';
    if (e.find( "Error" ) != std::string::npos) { err = true; underline = '^'; }
    else if (e.find( "Warning" ) != std::string::npos) underline = '~';

    // underline error and warning differently
    if (underline == '^' || underline == '~') {
      auto sloc = e.find( "at 1," );
      if (sloc != std::string::npos) {  // if we have location info
        // skip "at 1,"
        sloc += 5;
        // find a dot starting from after "at 1,"
        const auto eloc = e.find_first_of( '.', sloc );
        // compute location of error by extracting it from error message
        const std::size_t errloc = std::stoi( e.substr( sloc, eloc-sloc ) ) - 1;
        // find beginning of erroneous argument
        sloc = m_string.rfind( ' ', errloc-1 );
        // special-handle the first argument which has no space in front of it
        if (sloc == std::string::npos) sloc = 0; else ++sloc;
        std::cout << "s:" << sloc << ", e:" << errloc << std::endl;
        // underline error
        for (std::size_t i=sloc; i<errloc; ++i) pos[i] = underline;
      }
    }
  }

  // Output errors underlined (if any) to quiet stream, list errors, and exit
  if (!errors.empty()) {
    print % '\n';
    print % ">>> Parsed command line: '" % m_string % "'\n";
    print % ">>>                       " % pos % "\n";
    for (const auto& e : errors) print % ">>> " % e % std::endl;   // messages
    // Exit if there were any errors
    if (err) Throw( "Error(s) occurred while parsing the command line\n" );
  }
}
