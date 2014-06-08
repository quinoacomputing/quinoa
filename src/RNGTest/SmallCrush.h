//******************************************************************************
/*!
  \file      src/RNGTest/SmallCrush.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:15:19 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

#include <Option.h>
#include <Battery.h>
#include <TestU01Stack.h>

namespace rngtest {

//! SmallCrush
class SmallCrush {

  public:
    //! Return string identifying test suite name
    const std::string& name() const { return
      tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::SMALLCRUSH );
    }

    //! Add statistical tests to battery
    void addTests( std::vector< StatTest >& tests, tk::ctr::RNGType r );

  private:
    const TestU01Stack stack;
};

} // rngtest::

#endif // SmallCrush_h
