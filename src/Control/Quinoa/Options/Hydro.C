//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Hydro.C
  \author    J. Bakosi
  \date      Mon 13 Jan 2014 07:28:52 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's hydrodynamics model options
  \details   Quinoa's hydrodynamics model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Hydro.h>

using quinoa::ctr::Hydro;

// void
// Hydro::initFactory( HydroFactory& factory,
//                     std::list< HydroType >& reg ) const
// //******************************************************************************
// //  Register hydrodynamics models into factory
// //! \author  J. Bakosi
// //******************************************************************************
// {
//  reg.push_back( add< quinoa::SLM >( factory, HydroType::SLM ) );
//  reg.push_back( add< quinoa::GLM >( factory, HydroType::GLM ) );
// }
