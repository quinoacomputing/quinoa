//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Sun 06 Jul 2014 08:12:20 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver random number test suite driver
  \details   Driver random number test suite driver
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <Battery.h>
#include <RNGTestPrint.h>

//! Everything that contributes to the rngtest executable
namespace rngtest {

//! Battery factory type
using BatteryFactory = std::map< ctr::BatteryType, std::function< Battery() > >;

//! Random number generator test suite driver used polymorphically with Driver
class RNGTestDriver {

  public:
    //! Constructor
    explicit RNGTestDriver( int argc, char** argv, const RNGTestPrint& print );

    //! Execute driver
    void execute();

  private:
    const RNGTestPrint& m_print;
};

} // rngtest::

#endif // RNGTestDriver_h
