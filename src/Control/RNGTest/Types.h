//******************************************************************************
/*!
  \file      src/Control/RNGTest/Types.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 06:07:16 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Types for RNGTest's parsers
  \details   Types for RNGTest's parsers. This file defines the components of
    the tagged tuple that stores heteroegeneous objects in a hierarchical way.
    These components are therefore part of the grammar stack that is filled
    during parsing (both command-line argument parsing and control file
    parsing).
*/
//******************************************************************************
#ifndef RNGTestTypes_h
#define RNGTestTypes_h

#include <Tags.h>
#include <Types.h>
#include <RNGParam.h>
#include <RNGTest/Options/Battery.h>
#include <Options/RNG.h>

namespace rngtest {
namespace ctr {

using tk::ctr::RNGType;

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::battery, BatteryType,                   //!< Battery
  tag::rng,     std::vector< RNGType >         //!< Random number generators
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,   std::string                  //!< Control filename
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  #ifdef HAS_MKL
  tag::rngmkl,    tk::ctr::RNGMKLParameters,   //!< MKL RNG parameters
  #endif
  tag::rngsse,    tk::ctr::RNGSSEParameters    //!< RNGSSE parameters
>;

//! PEGTL location type to use throughout all of RNGTest's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // rngtest::

#endif // RNGTestTypes_h
