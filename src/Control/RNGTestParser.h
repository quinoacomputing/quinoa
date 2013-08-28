//******************************************************************************
/*!
  \file      src/Control/RNGTestParser.h
  \author    J. Bakosi
  \date      Wed Aug 28 15:18:05 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite parser
  \details   Random number generator test suite parser
*/
//******************************************************************************
#ifndef RNGTestParser_h
#define RNGTestParser_h

#include <fstream>
#include <vector>

#include <Parser.h>

namespace Quinoa {

//! RNGTestParser : Parser
class RNGTestParser : public Parser {

  public:
    //! Constructor
    //! \param[in]  filename  Control file name to read from
    //! \param[in]  control   Control object to put parsed data in
    explicit RNGTestParser(const std::string& filename, Control* const control)
      : Parser(filename, control) {}

    //! Destructor
    virtual ~RNGTestParser() noexcept = default;

    //! Parse random number generator test suite control file
    virtual void parse();

    //! Echo parsed information from random number generator test suite control
    virtual void echo() const;

  private:
    //! Don't permit copy constructor
    RNGTestParser(const RNGTestParser&) = delete;
    //! Don't permit copy assigment
    RNGTestParser& operator=(const RNGTestParser&) = delete;
    //! Don't permit move constructor
    RNGTestParser(RNGTestParser&&) = delete;
    //! Don't permit move assigment
    RNGTestParser& operator=(RNGTestParser&&) = delete;
};

} // namespace Quinoa

#endif // RNGTestParser_h
