//******************************************************************************
/*!
  \file      src/Control/RNGTestParser.h
  \author    J. Bakosi
  \date      Thu 26 Sep 2013 08:30:32 PM MDT
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

//! RNGTestParser : Parser
class RNGTestParser : public quinoa::Parser {

  public:
    //! Constructor
    //! \param[in]  filename  Control file name to read from
    //! \param[in]  control   Control object to put parsed data in
    explicit RNGTestParser(const std::string& filename,
                           RNGTestControl& control)
      : quinoa::Parser(filename),
        m_control(control) {}

    //! Destructor
    ~RNGTestParser() noexcept override = default;

    //! Parse random number generator test suite control file
    void parse() override;

  private:
    //! Don't permit copy constructor
    RNGTestParser(const RNGTestParser&) = delete;
    //! Don't permit copy assigment
    RNGTestParser& operator=(const RNGTestParser&) = delete;
    //! Don't permit move constructor
    RNGTestParser(RNGTestParser&&) = delete;
    //! Don't permit move assigment
    RNGTestParser& operator=(RNGTestParser&&) = delete;

    RNGTestControl& m_control;    //!< Control object
};

} // namespace rngtest

#endif // RNGTestParser_h
