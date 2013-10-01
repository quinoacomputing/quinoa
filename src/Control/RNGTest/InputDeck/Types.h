//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Types.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 10:19:19 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for random number generator test suite parsing
  \details   Types for random number generator test suite parsing
*/
//******************************************************************************
#ifndef RNGTestInputDeckTypes_h
#define RNGTestInputDeckTypes_h

#include <vector>
#include <string>
#include <tuple>

#include <RNGTest/Options/Battery.h>
#include <Quinoa/Options/RNG.h>

namespace rngtest {
//! control and parsing
namespace ctr {

//! Tags for Control's tagged tuple
struct title {};
struct suite {};
struct generator {};

} // ctr::
} // rngtest::

#endif // RNGTestInputDeckTypes_h
