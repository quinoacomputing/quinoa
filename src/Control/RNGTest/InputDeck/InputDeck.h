//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 03:28:47 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite input deck
  \details   Random number generator test suite input deck
*/
//******************************************************************************
#ifndef RNGTestInputDeck_h
#define RNGTestInputDeck_h

#include <Control.h>
#include <Option.h>
#include <Quinoa/Options/RNG.h>
#include <RNGTest/InputDeck/Types.h>

namespace rngtest {

//! InputDeck : Control< specialized to RNGTest >, see Types.h
//! This is also where the command line parser stores
class InputDeck :
  public quinoa::Control< // tag          type
                          ctr::title,     std::string,
                          ctr::selected,  ctr::selects,
                          ctr::io,        ctr::ios,
                          ctr::generator, std::vector<quinoa::sel::RNGType> > {

  public:
    //! Constructor: set defaults
    InputDeck() {
      using namespace ctr;
      // Default title
      set<title>("");
      // Default selections
      set<selected,battery>(sel::BatteryType::NO_BATTERY);
      // Default I/O parameters
      set<io,control>("");
      // Default requested generators
      set<generator>(std::vector<quinoa::sel::RNGType>());
    }

    //! Destructor
    ~InputDeck() noexcept override = default;

  private:
    //! Don't permit copy constructor
    InputDeck(const InputDeck&) = delete;
    //! Don't permit copy assigment
    InputDeck& operator=(const InputDeck&) = delete;
    //! Don't permit move constructor
    InputDeck(InputDeck&&) = delete;
    //! Don't permit move assigment
    InputDeck& operator=(InputDeck&&) = delete;
};

//! InputDeck defaults
static const InputDeck RNGTestDefaults;

} // namespace rngtest

#endif // RNGTestInputDeck_h
