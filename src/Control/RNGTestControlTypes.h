//******************************************************************************
/*!
  \file      src/Control/RNGTestControlTypes.h
  \author    J. Bakosi
  \date      Thu Sep 19 10:40:05 2013
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
//! control and parsing
namespace ctr {

//! Tags for Control's tagged tuple
struct title {};
struct suite {};
struct generator {};

} // ctr::
} // rngtest::

#endif // RNGTestControlTypes_h
