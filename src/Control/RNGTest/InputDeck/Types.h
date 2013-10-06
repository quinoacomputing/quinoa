//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Types.h
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 03:25:18 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for RNGTest's input deck parsing
  \details   Types for RNGTest's input deck parsing
*/
//******************************************************************************
#ifndef RNGTestInputDeckTypes_h
#define RNGTestInputDeckTypes_h

#include <RNGTest/InputDeck/Tags.h>
#include <RNGTest/Options/Battery.h>

namespace rngtest {
//! control and parsing
namespace ctr {

//! Storage of selected options
using selects = quinoa::tuple::tagged_tuple<
  battery,   sel::BatteryType     //!< Selected battery
>;

//! IO parameters storage
using ios = quinoa::tuple::tagged_tuple<
  control,     std::string       //!< Control filename
>;

} // ctr::
} // rngtest::

#endif // RNGTestInputDeckTypes_h
