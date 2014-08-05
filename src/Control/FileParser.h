//******************************************************************************
/*!
  \file      src/Control/FileParser.h
  \author    J. Bakosi
  \date      Mon 04 Aug 2014 07:23:23 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     File parser
  \details   File parser
*/
//******************************************************************************
#ifndef FileParser_h
#define FileParser_h

#include <string>

namespace tk {

//! FileParser
class FileParser {

  protected:
    //! Constructor
    explicit FileParser( const std::string& filename );

    const std::string m_filename;             //!< Name of file to parse
};

} // tk::

#endif // FileParser_h
