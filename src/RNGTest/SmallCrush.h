//******************************************************************************
/*!
  \file      src/RNGTest/SmallCrush.h
  \author    J. Bakosi
  \date      Sat 28 Jun 2014 02:52:18 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

#include <Option.h>
#include <Battery.h>
#include <StatTest.h>
#include <testu01suite.decl.h>

namespace rngtest {

//! SmallCrush
class SmallCrush {

  public:
    //! Return string identifying test suite name
    std::string name() const { return
      tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::SMALLCRUSH );
    }

    //! Add statistical tests to battery
    void addTests( std::vector< std::function< StatTest() > >& tests,
                   tk::ctr::RNGType rng,
                   CProxy_TestU01Suite& proxy );
};

} // rngtest::

#endif // SmallCrush_h
