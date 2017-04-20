// *****************************************************************************
/*!
  \file      src/Control/StringParser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     String parser base class declaration
  \details   String parser base class declaration. String parser base serves as
    a base class for various string parsers, e.g., command-line parsers. It does
    generic after-parser diagnostics.
*/
// *****************************************************************************
#ifndef StringParser_h
#define StringParser_h

#include <iosfwd>
#include <vector>

#include "Print.h"

namespace tk {

//! StringParser
class StringParser {

  protected:
    //! Constructor from C++-style std::string
    //! \param[in] string String to be parsed
    //! \author J. Bakosi
    explicit StringParser( const std::string& string ) : m_string( string ) {}

    //! Constructor from C-style command-line argument string
    explicit StringParser( int argc, char** argv );

    //! \brief Echo errors and warnings accumulated during parsing
    void diagnostics( const tk::Print& print,
                      const std::vector< std::string >& messages );

    std::string m_string;                     //!< String to parse
};

} // tk::

#endif // StringParser_h
