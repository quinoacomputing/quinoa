// *****************************************************************************
/*!
  \file      src/Control/FileParser.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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
                      const std::vector< std::string >& messages );

    const std::string m_filename;             //!< Name of file to parse
};

} // tk::

#endif // FileParser_h
