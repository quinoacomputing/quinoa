//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Hydro.C
  \author    J. Bakosi
  \date      Tue Oct 29 15:44:59 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's hydrodynamics model options
  \details   Quinoa's hydrodynamics model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Hydro.h>
#include <Hydro/SLM/SLM.h>
#include <Hydro/GLM/GLM.h>

using namespace quinoa::ctr;

void
Hydro::initFactory(HydroFactory& f, std::list<std::string>& reg) const
//******************************************************************************
//  Register hydrodynamics models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 reg.push_back( add<SLM>(f, HydroType::SLM) );
 reg.push_back( add<GLM>(f, HydroType::GLM) );
}
