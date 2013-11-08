//******************************************************************************
/*!
  \file      src/Control/RNGTest/Options/Battery.C
  \author    J. Bakosi
  \date      Thu Nov  7 11:34:57 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest's battery options
  \details   RNGTest's battery options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <RNGTest/Options/Battery.h>
#include <SmallCrush.h>
#include <Crush.h>
#include <BigCrush.h>

using namespace rngtest::ctr;

void
Battery::initFactory(BatteryFactory& f, std::list< BatteryType >& reg) const
//******************************************************************************
//  Register batteries into factory
//! \author  J. Bakosi
//******************************************************************************
{
  reg.push_back( add<SmallCrush>(f, BatteryType::SMALLCRUSH) );
  reg.push_back( add<Crush>(f, BatteryType::CRUSH) );
  reg.push_back( add<BigCrush>(f, BatteryType::BIGCRUSH) );
}
