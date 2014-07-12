//******************************************************************************
/*!
  \file      src/Control/StringParser.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:54:16 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     String parser
  \details   String parser
*/
//******************************************************************************
#ifndef StringParser_h
#define StringParser_h

#include <string>

#include <Parser.h>

namespace tk {

//! StringParser : Parser
class StringParser : public Parser {

  protected:
    //! Constructor from std::string
    explicit StringParser(const std::string& string) : m_string(string) {}

    //! Constructor from char**
    explicit StringParser(int argc, char** argv);

    std::string m_string;                     //!< String to parse
};

} // tk::

#endif // StringParser_h
