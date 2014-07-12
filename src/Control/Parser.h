//******************************************************************************
/*!
  \file      src/Control/Parser.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:54:54 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************
#ifndef Parser_h
#define Parser_h

namespace tk {

//! Parser base
class Parser {

  protected:
    //! Constructor
    explicit Parser() = default;

    //! Destructor
    virtual ~Parser() = default;
};

} // namespace tk

#endif // Parser_h
