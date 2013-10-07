//******************************************************************************
/*!
  \file      src/Control/StringParser.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:14:29 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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

    //! Destructor
    ~StringParser() noexcept override = default;

    //! StringParser's parse interface
    void parse() override = 0;

    std::string m_string;                     //!< String to parse

  private:
    //! Don't permit copy constructor
    StringParser(const StringParser&) = delete;
    //! Don't permit copy assigment
    StringParser& operator=(const StringParser&) = delete;
    //! Don't permit move constructor
    StringParser(StringParser&&) = delete;
    //! Don't permit move assigment
    StringParser& operator=(StringParser&&) = delete;
};

} // tk::

#endif // StringParser_h
