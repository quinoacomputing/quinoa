//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.C
  \author    J. Bakosi
  \date      Thu 21 Nov 2013 06:38:19 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************

#include <iostream>

extern "C" {
  #include <unif01.h>
  #include <bbattery.h>
  #include <swrite.h>
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
  ctr::RNGType r = m_base.control.get< ctr::selected, ctr::rng >()[0];
  if (r != ctr::RNGType::NO_RNG) {
    // Save name of random number generator
    quinoa::ctr::RNG rng;
    m_rngname = rng.name( r );
    // Instantiate random number generator
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
  swrite_Basic = FALSE; // no putput from TestU01

  unif01_Gen *gen;
  gen = unif01_CreateExternGen01( const_cast<char*>( m_rngname.c_str() ),
                                  MKLRNGUniform );
  bbattery_SmallCrush( gen );
  unif01_DeleteExternGen01( gen );
}
