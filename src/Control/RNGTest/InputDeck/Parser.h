//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Wed Oct  2 15:50:39 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite input deck parser
  \details   Random number generator test suite input deck parser
*/
//******************************************************************************
#ifndef RNGTestInputDeckParser_h
#define RNGTestInputDeckParser_h

#include <FileParser.h>
#include <Base.h>

namespace rngtest {

//! InputDeckParser : FileParser
class InputDeckParser : public quinoa::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser(const std::string& filename, quinoa::Base& base)
      : quinoa::FileParser(filename, base) {}

    //! Destructor
    ~InputDeckParser() noexcept override = default;

    //! Parse random number generator test suite control file
    void parse() override;

  private:
    //! Don't permit copy constructor
    InputDeckParser(const InputDeckParser&) = delete;
    //! Don't permit copy assigment
    InputDeckParser& operator=(const InputDeckParser&) = delete;
    //! Don't permit move constructor
    InputDeckParser(InputDeckParser&&) = delete;
    //! Don't permit move assigment
    InputDeckParser& operator=(InputDeckParser&&) = delete;
};

} // namespace rngtest

#endif // RNGTestInputDeckParser_h
