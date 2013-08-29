//******************************************************************************
/*!
  \file      src/Control/Parser.h
  \author    J. Bakosi
  \date      Wed 28 Aug 2013 08:46:32 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************
#ifndef Parser_h
#define Parser_h

#include <string>

namespace Quinoa {

//! Parser base
class Parser {

  public:
    //! Constructor
    explicit Parser(const std::string& filename);

    //! Destructor
    virtual ~Parser() noexcept = default;

    //! Parser interface
    virtual void parse() = 0;

    //! Echo parsed information interface
    virtual void echo() const = 0;

  protected:
    const std::string m_filename;                     //!< Name of file to parse

  private:
    //! Don't permit copy constructor
    Parser(const Parser&) = delete;
    //! Don't permit copy assigment
    Parser& operator=(const Parser&) = delete;
    //! Don't permit move constructor
    Parser(Parser&&) = delete;
    //! Don't permit move assigment
    Parser& operator=(Parser&&) = delete;
};

} // namespace Quinoa

#endif // Parser_h
