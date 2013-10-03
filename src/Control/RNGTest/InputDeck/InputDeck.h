//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Thu Oct  3 17:27:32 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite input deck
  \details   Random number generator test suite input deck
*/
//******************************************************************************
#ifndef RNGTestInputDeck_h
#define RNGTestInputDeck_h

#include <Control.h>
#include <RNGTest/InputDeck/Types.h>

namespace rngtest {

//! InputDeck : Control< specialized to RNGTest >, see Types.h
//! This is also where the command line parser stores
class InputDeck :
  public quinoa::Control< // tag          type
                          ctr::title,     std::string,
                          ctr::suite,     sel::BatteryType,
                          ctr::generator, std::vector<quinoa::sel::RNGType> > {

  public:
    //! Constructor: set defaults
    explicit InputDeck() = default;

// //! Default bundle for RNGTest's control
// const Bundle defaults(
//   "",                                  //!< Title
//   sel::BatteryType::NO_BATTERY,     //!< RNG test suite
//   std::vector<quinoa::sel::RNGType>()       //!< Random number generators
// );

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

} // namespace rngtest

#endif // RNGTestInputDeck_h
