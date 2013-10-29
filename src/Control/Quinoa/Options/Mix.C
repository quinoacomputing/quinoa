//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mix.C
  \author    J. Bakosi
  \date      Tue Oct 29 15:45:18 2013
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
Mix::initFactory(MixFactory& f, std::list<std::string>& reg) const
//******************************************************************************
//  Register material mix models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 reg.push_back( add<Dirichlet>(f, MixType::DIRICHLET) );
 reg.push_back( add<GenDirichlet>(f, MixType::GENERALIZED_DIRICHLET) );
}
