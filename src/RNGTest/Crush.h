//******************************************************************************
/*!
  \file      src/RNGTest/Crush.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:15:34 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Crush battery
  \details   Crush battery
*/
//******************************************************************************
#ifndef Crush_h
#define Crush_h

#include <Option.h>
#include <Battery.h>
#include <TestU01Stack.h>

namespace rngtest {

//! Crush : TestU01Suite
class Crush {

  public:
    //! Return string identifying test suite name
    const std::string& name() const { return
      tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::CRUSH );
    }

    //! Add statistical tests to battery
    void addTests( std::vector< StatTest >& tests, tk::ctr::RNGType r );

  private:
    const TestU01Stack stack;
};

} // rngtest::

#endif // Crush_h
