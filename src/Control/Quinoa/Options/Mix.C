//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mix.C
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:46:51 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's material mix model options
  \details   Quinoa's material mix model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Mix.h>
#include <DM.h>
#include <GDM.h>

using quinoa::ctr::Mix;

void
Mix::initFactory( MixFactory& factory, std::list< MixType >& reg ) const
//******************************************************************************
//  Register material mix models into factory
//! \author  J. Bakosi
//******************************************************************************
{
 reg.push_back( add< DM >( factory, MixType::DM ) );
 reg.push_back( add< GDM >( factory, MixType::GDM ) );
}
