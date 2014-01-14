//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Mix.C
  \author    J. Bakosi
  \date      Mon 13 Jan 2014 07:29:19 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's material mix model options
  \details   Quinoa's material mix model options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/Mix.h>

using quinoa::ctr::Mix;

// void
// Mix::initFactory( MixFactory& factory, std::list< MixType >& reg ) const
// //******************************************************************************
// //  Register material mix models into factory
// //! \author  J. Bakosi
// //******************************************************************************
// {
//  reg.push_back( add< DM >( factory, MixType::DM ) );
//  reg.push_back( add< GDM >( factory, MixType::GDM ) );
// }
