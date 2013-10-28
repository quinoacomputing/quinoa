//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Hydro.C
  \author    J. Bakosi
  \date      Mon Oct 28 08:52:15 2013
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
Hydro::initFactory(HydroFactory& f) const
//******************************************************************************
//  Register hydrodynamics models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 f[ HydroType::SLM ] = boost::factory< SLM* >();
 f[ HydroType::GLM ] = boost::factory< GLM* >();
}
