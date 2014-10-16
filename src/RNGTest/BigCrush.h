//******************************************************************************
/*!
  \file      src/RNGTest/BigCrush.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 09:50:21 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************
#ifndef BigCrush_h
#define BigCrush_h

#include <Battery.h>
#include <StatTest.h>
#include <testu01suite.decl.h>

namespace rngtest {

//! BigCrush
class BigCrush {

  public:
    //! Return string identifying test suite name
    std::string name() const
    { return ctr::Battery().name( rngtest::ctr::BatteryType::BIGCRUSH ); }

    //! Add statistical tests to battery
    void addTests( std::vector< std::function< StatTest() > >& tests,
                   tk::ctr::RNGType rng,
                   CProxy_TestU01Suite& proxy );
};

} // rngtest::

#endif // BigCrush_h
