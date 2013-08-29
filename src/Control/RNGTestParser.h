//******************************************************************************
/*!
  \file      src/Control/RNGTestParser.h
  \author    J. Bakosi
  \date      Thu Aug 29 16:58:22 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite parser
  \details   Random number generator test suite parser
*/
//******************************************************************************
#ifndef RNGTestParser_h
#define RNGTestParser_h

#include <Parser.h>
#include <RNGTestControl.h>

namespace rngtest {

using quinoa::Parser;

//! RNGTestParser : Parser
class RNGTestParser : public Parser {

  public:
    //! Constructor
    //! \param[in]  filename  Control file name to read from
    //! \param[in]  control   Control object to put parsed data in
    explicit RNGTestParser(const std::string& filename,
                           RNGTestControl* const control)
      : Parser(filename),
        m_control(control) {}

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

    RNGTestControl* const m_control;    //!< Control object
};

} // namespace rngtest

#endif // RNGTestParser_h
