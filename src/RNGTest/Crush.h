//******************************************************************************
/*!
  \file      src/RNGTest/Crush.h
  \author    J. Bakosi
  \date      Sat 28 Jun 2014 02:53:45 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Crush battery
  \details   Crush battery
*/
//******************************************************************************
#ifndef Crush_h
#define Crush_h

#include <Option.h>
#include <Battery.h>
#include <StatTest.h>
#include <testu01suite.decl.h>

namespace rngtest {

//! Crush
class Crush {

  public:
    //! Return string identifying test suite name
    std::string name() const { return
      tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::CRUSH );
    }

    //! Add statistical tests to battery
    void addTests( std::vector< std::function< StatTest() > >& tests,
                   tk::ctr::RNGType rng,
                   CProxy_TestU01Suite& proxy );
};

} // rngtest::

#endif // Crush_h
