//******************************************************************************
/*!
  \file      src/RNGTest/Crush.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 04:24:03 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Class re-creating the TestU01 library's Crush battery
  \details   Class re-creating the TestU01 library's Crush battery.
*/
//******************************************************************************
#ifndef Crush_h
#define Crush_h

#include <Battery.h>
#include <StatTest.h>
#include <testu01suite.decl.h>

namespace rngtest {

//! Class registering the TestU01 library's Crush battery
class Crush {

  public:
    //! Return string identifying test suite name
    //! \return Test suite name as a std::string
    std::string name() const
    { return ctr::Battery().name( rngtest::ctr::BatteryType::CRUSH ); }

    //! Add statistical tests to battery
    void addTests( std::vector< std::function< StatTest() > >& tests,
                   tk::ctr::RNGType rng,
                   CProxy_TestU01Suite& proxy );
};

} // rngtest::

#endif // Crush_h
