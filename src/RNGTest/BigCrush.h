//******************************************************************************
/*!
  \file      src/RNGTest/BigCrush.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:15:48 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************
#ifndef BigCrush_h
#define BigCrush_h

#include <Option.h>
#include <Battery.h>
#include <TestU01Stack.h>

namespace rngtest {

//! BigCrush
class BigCrush {

  public:
    //! Return string identifying test suite name
    const std::string& name() const { return
      tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::BIGCRUSH );
    }

    //! Add statistical tests to battery
    void addTests( std::vector< StatTest >& tests, tk::ctr::RNGType r );

  private:
    const TestU01Stack stack;
};

} // rngtest::

#endif // BigCrush_h
