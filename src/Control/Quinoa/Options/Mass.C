//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mass.C
  \author    J. Bakosi
  \date      Mon Oct 28 08:47:17 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's mass model options
  \details   Quinoa's mass model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Mass.h>
#include <Mass/Beta/Beta.h>

using namespace quinoa::ctr;


void
Mass::initFactory(MassFactory& f) const
//******************************************************************************
//  Register mass models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 f[ MassType::BETA ] = boost::factory< Beta* >();
}
