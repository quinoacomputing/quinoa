//******************************************************************************
/*!
  \file      src/Control/FileParser.h
  \author    J. Bakosi
  \date      Thu 28 Aug 2014 03:27:45 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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
    void diagnostics( const tk::Print& print,
                      const std::vector< std::string >& errors );

    const std::string m_filename;             //!< Name of file to parse
};

} // tk::

#endif // FileParser_h
