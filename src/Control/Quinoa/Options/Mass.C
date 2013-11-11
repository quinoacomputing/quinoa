//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mass.C
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 10:33:54 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's mass model options
  \details   Quinoa's mass model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Mass.h>
#include <Mass/Beta/Beta.h>

using quinoa::ctr::Mass;

void
Mass::initFactory( MassFactory& factory, std::list< MassType >& reg ) const
//******************************************************************************
//  Register mass models into factory
//! \author  J. Bakosi
//******************************************************************************
{
  reg.push_back( add< Beta >( factory, MassType::BETA ) );
}
