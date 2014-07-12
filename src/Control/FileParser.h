//******************************************************************************
/*!
  \file      src/Control/FileParser.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:53:57 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     File parser
  \details   File parser
*/
//******************************************************************************
#ifndef FileParser_h
#define FileParser_h

#include <string>

#include <Parser.h>

namespace tk {

//! FileParser : Parser
class FileParser : public Parser {

  protected:
    //! Constructor
    explicit FileParser(const std::string& filename);

    const std::string m_filename;             //!< Name of file to parse
};

} // tk::

#endif // FileParser_h
