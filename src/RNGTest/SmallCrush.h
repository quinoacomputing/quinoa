//******************************************************************************
/*!
  \file      src/RNGTest/SmallCrush.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 09:49:33 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

#include <Battery.h>
#include <StatTest.h>
#include <testu01suite.decl.h>

namespace rngtest {

//! SmallCrush
class SmallCrush {

  public:
    //! Return string identifying test suite name
    std::string name() const
    { return ctr::Battery().name( rngtest::ctr::BatteryType::SMALLCRUSH ); }

    //! Add statistical tests to battery
    void addTests( std::vector< std::function< StatTest() > >& tests,
                   tk::ctr::RNGType rng,
                   CProxy_TestU01Suite& proxy );
};

} // rngtest::

#endif // SmallCrush_h
