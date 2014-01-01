//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mass.C
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:43:26 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's mass model options
  \details   Quinoa's mass model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Mass.h>
#include <BM.h>

using quinoa::ctr::Mass;

void
Mass::initFactory( MassFactory& factory, std::list< MassType >& reg ) const
//******************************************************************************
//  Register mass models into factory
//! \author  J. Bakosi
//******************************************************************************
{
  reg.push_back( add< BM >( factory, MassType::BM ) );
}
