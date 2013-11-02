//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mass.C
  \author    J. Bakosi
  \date      Sat 02 Nov 2013 10:51:23 AM MDT
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
Mass::initFactory(MassFactory& f, std::list<std::string>& reg) const
//******************************************************************************
//  Register mass models into factory
//! \author  J. Bakosi
//******************************************************************************
{
  reg.push_back( add<Beta>(f, MassType::BETA) );
}
