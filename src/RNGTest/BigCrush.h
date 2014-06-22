//******************************************************************************
/*!
  \file      src/RNGTest/BigCrush.h
  \author    J. Bakosi
  \date      Mon 16 Jun 2014 11:22:12 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************
#ifndef BigCrush_h
#define BigCrush_h

#include <Option.h>
#include <Battery.h>
#include <StatTest.h>
#include <testu01suite.decl.h>

namespace rngtest {

//! BigCrush
class BigCrush {

  public:
    //! Return string identifying test suite name
    const std::string& name() const { return
      tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::BIGCRUSH );
    }

    //! Add statistical tests to battery
    void addTests( std::vector< std::function< StatTest() > >& tests,
                   tk::ctr::RNGType rng,
                   CProxy_TestU01Suite& proxy );
};

} // rngtest::

#endif // BigCrush_h
