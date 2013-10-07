//******************************************************************************
/*!
  \file      src/Control/Parser.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:14:11 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************
#ifndef Parser_h
#define Parser_h

#include <Base.h>

namespace tk {

//! Parser base
class Parser {

  protected:
    //! Constructor
    explicit Parser() = default;

    //! Destructor
    virtual ~Parser() noexcept = default;

    //! Parse interface
    virtual void parse() = 0;

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

} // namespace tk

#endif // Parser_h
