//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mass.C
  \author    J. Bakosi
  \date      Thu Nov  7 11:33:13 2013
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
Mass::initFactory(MassFactory& f, std::list< MassType >& reg) const
//******************************************************************************
//  Register mass models into factory
//! \author  J. Bakosi
//******************************************************************************
{
  reg.push_back( add<Beta>(f, MassType::BETA) );
}
