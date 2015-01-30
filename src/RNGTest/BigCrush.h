//******************************************************************************
/*!
  \file      src/RNGTest/BigCrush.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 04:22:30 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Class re-creating the TestU01 library's BigCrush battery
  \details   Class re-creating the TestU01 library's BigCrush battery.
*/
//******************************************************************************
#ifndef BigCrush_h
#define BigCrush_h

#include <Battery.h>
#include <StatTest.h>
#include <testu01suite.decl.h>

namespace rngtest {

//! Class registering the TestU01 library's BigCrush battery
class BigCrush {

  public:
    //! Return string identifying test suite name
    //! \return Test suite name as a std::string
    std::string name() const
    { return ctr::Battery().name( rngtest::ctr::BatteryType::BIGCRUSH ); }

    //! Add statistical tests to battery
    void addTests( std::vector< std::function< StatTest() > >& tests,
                   tk::ctr::RNGType rng,
                   CProxy_TestU01Suite& proxy );
};

} // rngtest::

#endif // BigCrush_h
