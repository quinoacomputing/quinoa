//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Sun 03 Nov 2013 08:46:08 PM MST
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
  public tk::Control< // tag      type
                      title,      std::string,
                      selected,   selects,
                      io,         ios,
                      cmd,        CmdLine,
                      param,      parameters > {

  public:
    //! Constructor: set defaults
    InputDeck() {
      using namespace ctr;
      // Default title
      set< title >( "" );
      // Default selections
      set< selected, battery >( BatteryType::NO_BATTERY );
      set< selected, mkl_uniform_method >( MKLUniformMethodType::NO_METHOD );
      // Default I/O parameters
      set< io, control >( "" );
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
