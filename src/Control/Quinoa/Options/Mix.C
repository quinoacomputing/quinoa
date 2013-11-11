//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mix.C
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 10:35:31 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's material mix model options
  \details   Quinoa's material mix model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Mix.h>
#include <Mix/Dirichlet/Dirichlet.h>
#include <Mix/GenDirichlet/GenDirichlet.h>

using quinoa::ctr::Mix;

void
Mix::initFactory( MixFactory& factory, std::list< MixType >& reg ) const
//******************************************************************************
//  Register material mix models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 reg.push_back( add< Dirichlet >( factory, MixType::DIRICHLET ) );
 reg.push_back( add< GenDirichlet >( factory, MixType::GENERALIZED_DIRICHLET ) );
}
