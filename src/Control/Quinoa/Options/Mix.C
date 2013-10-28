//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mix.C
  \author    J. Bakosi
  \date      Mon Oct 28 08:51:03 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's material mix model options
  \details   Quinoa's material mix model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Mix.h>
#include <Mix/Dirichlet/Dirichlet.h>
#include <Mix/GenDirichlet/GenDirichlet.h>

using namespace quinoa::ctr;

void
Mix::initFactory(MixFactory& f) const
//******************************************************************************
//  Register material mix models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 f[ MixType::DIRICHLET ] = boost::factory< Dirichlet* >();
 f[ MixType::GENERALIZED_DIRICHLET ] = boost::factory< GenDirichlet* >();
}
