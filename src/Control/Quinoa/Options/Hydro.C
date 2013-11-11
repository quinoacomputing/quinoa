//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Hydro.C
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 10:31:43 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's hydrodynamics model options
  \details   Quinoa's hydrodynamics model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Hydro.h>
#include <Hydro/SLM/SLM.h>
#include <Hydro/GLM/GLM.h>

using quinoa::ctr::Hydro;

void
Hydro::initFactory( HydroFactory& factory,
                    std::list< HydroType >& reg ) const
//******************************************************************************
//  Register hydrodynamics models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 reg.push_back( add< quinoa::SLM >( factory, HydroType::SLM ) );
 reg.push_back( add< quinoa::GLM >( factory, HydroType::GLM ) );
}
