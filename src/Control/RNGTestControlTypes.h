//******************************************************************************
/*!
  \file      src/Control/RNGTestControlTypes.h
  \author    J. Bakosi
  \date      Tue 03 Sep 2013 10:48:26 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for random number generator test suite control and parsing
  \details   Types for random number generator test suite control and parsing
*/
//******************************************************************************
#ifndef RNGTestControlTypes_h
#define RNGTestControlTypes_h

#include <vector>
#include <string>
#include <tuple>

#include <RNGTestOptions.h>
#include <RNGOptions.h>

namespace rngtest {

namespace control {

//! Tags for Control's tagged tuple
struct title {};
struct suite {};
struct generator {};

} // namespace control

} // namespace rngtest

#endif // RNGTestControlTypes_h
