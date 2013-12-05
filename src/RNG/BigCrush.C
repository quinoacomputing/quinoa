//******************************************************************************
/*!
  \file      src/RNG/BigCrush.C
  \author    J. Bakosi
  \date      Wed 04 Dec 2013 09:17:47 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************

#include <BigCrush.h>
#include <TestU01.h>

using rngtest::BigCrush;

BigCrush::BigCrush(const Base& base) :
  TestU01Suite( base,
    tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::BIGCRUSH ) )
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  setupRNGs( *this );
}

void
BigCrush::addTests( const quinoa::ctr::RNGType& rng, const Gen01Ptr& gen )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
}
