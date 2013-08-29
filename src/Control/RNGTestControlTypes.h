//******************************************************************************
/*!
  \file      src/Control/RNGTestControlTypes.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:03:48 2013
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

//! Position enum for accessing fields of tuple Bundle using names as in struct
enum BundlePosition { TITLE=0,
                      RNGTEST,
                      RNGS
};

//! Storage bundle for parsed data
using Bundle = std::tuple<
  std::string,                  //!< Test suite title
  select::RNGTestType,          //!< Selected RNG test suite
  std::vector<select::RNGType>  //!< Random number generators
>;

//! Default bundle for RNGTest's control
const Bundle defaults(
  "",                                  //!< Title
  select::RNGTestType::NO_RNGTEST,     //!< RNG test suite
  std::vector<select::RNGType>()       //!< Random number generators
);


//! Vector of bools indicating whether data is set in Bundle during parsing
using BoolBundle = std::vector<bool>;

} // namespace control

} // namespace rngtest

#endif // RNGTestControlTypes_h
