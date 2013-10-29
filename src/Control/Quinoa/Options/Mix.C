//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mix.C
  \author    J. Bakosi
  \date      Tue Oct 29 15:36:57 2013
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
Mix::initFactory(MixFactory& f, std::list<std::string>& names) const
//******************************************************************************
//  Register material mix models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 names.push_back( add<Dirichlet>(f, MixType::DIRICHLET) );
 names.push_back( add<GenDirichlet>(f, MixType::GENERALIZED_DIRICHLET) );
}
