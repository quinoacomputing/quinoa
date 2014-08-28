//******************************************************************************
/*!
  \file      src/Control/StringParser.h
  \author    J. Bakosi
  \date      Thu 28 Aug 2014 12:53:36 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     String parser
  \details   String parser
*/
//******************************************************************************
#ifndef StringParser_h
#define StringParser_h

#include <string>

#include <Print.h>

namespace tk {

//! StringParser
class StringParser {

  protected:
    //! Constructor from std::string
    explicit StringParser( const std::string& string ) : m_string( string ) {}

    //! Constructor from char**
    explicit StringParser( int argc, char** argv );

    //! Echo errors and warnings accumulated during parsing
    void echoErrors( const tk::Print& print,
                     const std::vector< std::string >& messages );

    std::string m_string;                     //!< String to parse
};

} // tk::

#endif // StringParser_h
