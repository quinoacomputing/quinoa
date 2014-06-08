//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 01:36:49 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver random number test suite driver
  \details   Driver random number test suite driver
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <Battery.h>

//! Everything that contributes to the rngtest executable
namespace rngtest {

//! Battery factory type
using BatteryFactory = std::map< ctr::BatteryType, std::function< Battery() > >;

//! Random number generator test suite driver used polymorphically with Driver
class RNGTestDriver {

  public:
    //! Constructor
    explicit RNGTestDriver( int argc, char** argv );

    //! Execute driver
    void execute();

  private:
    ctr::InputDeck m_control;
};

} // rngtest::

#endif // RNGTestDriver_h
