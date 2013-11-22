//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.C
  \author    J. Bakosi
  \date      Thu 21 Nov 2013 06:16:00 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************

#include <iostream>

extern "C" {
  #include <unif01.h>
  #include <bbattery.h>
}

#include <SmallCrush.h>
#include <MKLRNGWrappers.h>

namespace rngtest {

extern std::unique_ptr< tk::RNG > g_rng;

} // rngtest::

using rngtest::SmallCrush;

SmallCrush::SmallCrush(const Base& base) : TestU01(base)
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate random number generator
  ctr::RNGType r = m_base.control.get< ctr::selected, ctr::rng >()[0];
  if (r != ctr::RNGType::NO_RNG) {
    g_rng = std::unique_ptr< tk::RNG >( m_base.rng[r]() );
  }
}

void
SmallCrush::run()
//******************************************************************************
//  Run SmallCrush battery
//! \author  J. Bakosi
//******************************************************************************
{
  const char* name = "RNG test name";

  unif01_Gen *gen;
  gen = unif01_CreateExternGen01( const_cast<char*>(name), MKLRNGUniform );
  bbattery_SmallCrush( gen );
  unif01_DeleteExternGen01( gen );
}
