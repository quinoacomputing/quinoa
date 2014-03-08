//******************************************************************************
/*!
  \file      src/Control/RNGTest/Types.h
  \author    J. Bakosi
  \date      Sat 08 Mar 2014 06:43:05 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for RNGTest's parsers
  \details   Types for RNGTest's parsers
*/
//******************************************************************************
#ifndef RNGTestTypes_h
#define RNGTestTypes_h

#include <RNGTest/Tags.h>
#include <RNGTest/Options/Battery.h>
#include <Options/RNG.h>

namespace rngtest {
//! control and parsing
namespace ctr {

using tk::ctr::RNGType;

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::battery,   BatteryType,                //!< Battery
  tk::tag::rng,   std::vector< RNGType >      //!< Random number generators
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,   std::string                 //!< Control filename
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  #ifdef HAS_MKL
  tk::tag::mklrng,    tk::ctr::MKLRNGParameters,   //!< MKL RNG parameters
  #endif
  tk::tag::rngsse,    tk::ctr::RNGSSEParameters    //!< RNGSSE parameters
>;

//! PEGTL location type to use throughout all of RNGTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // rngtest::

#endif // RNGTestTypes_h
