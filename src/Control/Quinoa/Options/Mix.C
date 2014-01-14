//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mix.C
  \author    J. Bakosi
  \date      Tue Jan 14 07:57:26 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's material mix model options
  \details   Quinoa's material mix model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Mix.h>

using quinoa::ctr::Mix;

void
Mix::initFactory( MixFactory& factory, std::list< MixType >& reg ) const
//******************************************************************************
//  Register material mix models into factory
//! \author  J. Bakosi
//******************************************************************************
{
// reg.push_back( add< Dirichlet< InitRaw, DirCoeffConst > >
//                   ( factory, MixType::DIRICHLET ) );
// reg.push_back( add< GDM >( factory, MixType::GDM ) );
}
