//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.C
  \author    J. Bakosi
  \date      Thu 21 Nov 2013 02:40:57 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************

#include <iostream>

#include <SmallCrush.h>

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
    m_rng = std::unique_ptr< tk::RNG >( m_base.rng[r]() );
  }
}

void
SmallCrush::run()
//******************************************************************************
//  Run SmallCrush battery
//! \author  J. Bakosi
//******************************************************************************
{
  
}
