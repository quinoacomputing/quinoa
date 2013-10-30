//******************************************************************************
/*!
  \file      src/Control/RNGTest/Types.h
  \author    J. Bakosi
  \date      Wed Oct 30 07:14:17 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for RNGTest's parsers
  \details   Types for RNGTest's parsers
*/
//******************************************************************************
#ifndef RNGTestTypes_h
#define RNGTestTypes_h

#include <RNGTest/Tags.h>
#include <RNGTest/Options/Battery.h>
#include <Quinoa/Options/RNG.h>

namespace rngtest {
//! control and parsing
namespace ctr {

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  battery,   ctr::BatteryType     //!< Selected battery
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  control,     std::string       //!< Control filename
>;

//! PEGTL location type to use throughout all of RNGTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // rngtest::

#endif // RNGTestTypes_h
