//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Wed 09 Oct 2013 10:28:28 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite input deck parser
  \details   Random number generator test suite input deck parser
*/
//******************************************************************************
#ifndef RNGTestInputDeckParser_h
#define RNGTestInputDeckParser_h

#include <FileParser.h>

namespace rngtest {

//! InputDeckParser : FileParser
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser(Base& base);

    //! Destructor
    ~InputDeckParser() noexcept override = default;

  private:
    //! Don't permit copy constructor
    InputDeckParser(const InputDeckParser&) = delete;
    //! Don't permit copy assigment
    InputDeckParser& operator=(const InputDeckParser&) = delete;
    //! Don't permit move constructor
    InputDeckParser(InputDeckParser&&) = delete;
    //! Don't permit move assigment
    InputDeckParser& operator=(InputDeckParser&&) = delete;

    const Base& m_base;                  //!< Essentials
};

} // namespace rngtest

#endif // RNGTestInputDeckParser_h
