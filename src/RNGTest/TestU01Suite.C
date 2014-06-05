//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.C
  \author    J. Bakosi
  \date      Wed 21 May 2014 03:17:06 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 suite
  \details   TestU01 suite
*/
//******************************************************************************

#include <TestU01Suite.h>
#include <SmallCrush.h>
#include <Crush.h>
#include <BigCrush.h>

#include <pup.h>

using rngtest::TestU01Suite;

//! Specialized TestU01Suite PUP ids

template<> PUP::able::PUP_ID TestU01Suite< rngtest::SmallCrush >::my_PUP_ID = 0;
template<> PUP::able::PUP_ID TestU01Suite< rngtest::Crush >::my_PUP_ID = 0;
template<> PUP::able::PUP_ID TestU01Suite< rngtest::BigCrush >::my_PUP_ID = 0;

#include <testu01suite.def.h>
