//******************************************************************************
/*!
  \file      src/Control/FileParser.h
  \author    J. Bakosi
  \date      Tue 26 Aug 2014 06:12:13 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     File parser
  \details   File parser
*/
//******************************************************************************
#ifndef FileParser_h
#define FileParser_h

#include <string>

#include <Print.h>

namespace tk {

//! FileParser
class FileParser {

  protected:
    //! Constructor
    explicit FileParser( const std::string& filename );

    //! Echo errors accumulated during parsing
    void echoErrors( const tk::Print& print,
                     const std::vector< std::string >& errors );

    const std::string m_filename;             //!< Name of file to parse
};

} // tk::

#endif // FileParser_h
