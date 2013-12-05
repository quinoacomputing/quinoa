//******************************************************************************
/*!
  \file      src/RNG/Crush.C
  \author    J. Bakosi
  \date      Wed 04 Dec 2013 09:23:16 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Crush battery
  \details   Crush battery
*/
//******************************************************************************

#include <Crush.h>
#include <TestU01.h>

using rngtest::Crush;

Crush::Crush(const Base& base) :
  TestU01Suite( base,
    tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::CRUSH ) )
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  setupRNGs( *this );
}

void
Crush::addTests( const quinoa::ctr::RNGType& rng, const Gen01Ptr& gen )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
}
