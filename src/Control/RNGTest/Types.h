//******************************************************************************
/*!
  \file      src/Control/RNGTest/Types.h
  \author    J. Bakosi
  \date      Sun 03 Nov 2013 08:26:54 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for RNGTest's parsers
  \details   Types for RNGTest's parsers
*/
//******************************************************************************
#ifndef RNGTestTypes_h
#define RNGTestTypes_h

#include <RNGTest/Tags.h>
#include <RNGTest/Options/Battery.h>
#include <RNGTest/Options/MKLUniformMethod.h>
#include <Quinoa/Options/RNG.h>

namespace rngtest {
//! control and parsing
namespace ctr {

using quinoa::ctr::RNGType;

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  battery,             ctr::BatteryType,           //!< Battery
  rng,                 std::vector< RNGType >,     //!< Random number generators
  mkl_uniform_method,  ctr::MKLUniformMethodType   //!< MKL uniform RNG method
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  control,   std::string                           //!< Control filename
>;

//! Random number generator parameters storage
using RNGParameters = tk::tuple::tagged_tuple<
  seed,                std::vector< unsigned int > //!< RNG seeds
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  rng,       RNGParameters                         //!< Random number generators
>;

//! PEGTL location type to use throughout all of RNGTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // rngtest::

#endif // RNGTestTypes_h
