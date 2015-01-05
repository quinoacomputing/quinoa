//******************************************************************************
/*!
  \file      src/Control/FileParser.h
  \author    J. Bakosi
  \date      Wed 14 Jan 2015 02:30:26 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     File parser base class declaration
  \details   File parser base class declaration. File parser base serves as a
    base class for various file parsers, e.g., input deck parsers. It does
    generic low-level I/O, e.g., testing whether the file to be parsed exits or
    not and associated error handling, as well as after-parser diagnostics.
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
    //! \brief Constructor
    explicit FileParser( const std::string& filename );

    //! \brief Echo errors accumulated during parsing
    void diagnostics( const tk::Print& print,
                      const std::vector< std::string >& errors );

    const std::string m_filename;             //!< Name of file to parse
};

} // tk::

#endif // FileParser_h
