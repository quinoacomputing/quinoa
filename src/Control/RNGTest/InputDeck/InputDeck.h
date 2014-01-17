//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 09:56:27 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite input deck
  \details   Random number generator test suite input deck
*/
//******************************************************************************
#ifndef RNGTestInputDeck_h
#define RNGTestInputDeck_h

#include <Control.h>
#include <Option.h>
#include <RNGTest/CmdLine/CmdLine.h>

namespace rngtest {
namespace ctr {

//! InputDeck : Control< specialized to RNGTest >, see Types.h
//! This is also where the command line parser stores
class InputDeck :
  public tk::Control< // tag           type
                      tag::title,      std::string,
                      tag::selected,   selects,
                      tag::io,         ios,
                      tag::cmd,        CmdLine,
                      tag::param,      parameters > {

  public:
    //! Constructor: set defaults
    InputDeck() {
      // Default title
      set< tag::title >( "" );
      // Default I/O parameters
      set< tag::io, tag::control >( "" );
    }

    //! Destructor
    ~InputDeck() noexcept override = default;

    //! Instruct compiler to generate copy assigment
    InputDeck& operator=(const InputDeck&) = default;

  private:
    //! Don't permit copy constructor
    InputDeck(const InputDeck&) = delete;
    //! Don't permit move constructor
    InputDeck(InputDeck&&) = delete;
    //! Don't permit move assigment
    InputDeck& operator=(InputDeck&&) = delete;
};

//! InputDeck defaults
static const InputDeck InputDeckDefaults;

} // ctr::
} // rngtest::

#endif // RNGTestInputDeck_h
