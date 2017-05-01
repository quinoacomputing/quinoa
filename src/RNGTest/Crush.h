// *****************************************************************************
/*!
  \file      src/RNGTest/Crush.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Class re-creating the TestU01 library's Crush battery
  \details   Class re-creating the TestU01 library's Crush battery.
*/
// *****************************************************************************
#ifndef Crush_h
#define Crush_h

#include <vector>
#include <functional>
#include <iosfwd>

#include "Options/RNG.h"
#include "RNGTest/Options/Battery.h"

namespace rngtest {

class CProxy_TestU01Suite;
class StatTest;

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
