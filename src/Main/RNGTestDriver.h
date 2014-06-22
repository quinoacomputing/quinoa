//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Mon 16 Jun 2014 07:26:46 PM MDT
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
};

} // rngtest::

#endif // RNGTestDriver_h
