// *****************************************************************************
/*!
  \file      src/Control/FileParser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     File parser base class declaration
  \details   File parser base class declaration. File parser base serves as a
    base class for various file parsers, e.g., input deck parsers. It does
    generic low-level I/O, e.g., testing whether the file to be parsed exits or
    not and associated error handling, as well as after-parser diagnostics.
*/
// *****************************************************************************
#ifndef FileParser_h
#define FileParser_h

#include <vector>
#include <string>

namespace tk {

class Print;

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
