// *****************************************************************************
/*!
  \file      tests/unit/Base/TestFlip_map.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/Flip_map.hpp
  \details   Unit tests for Base/Flip_map.hpp
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Flip_map.hpp"
#include "Types.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct Flip_map_common {};

//! Test group shortcuts
using Flip_map_group = test_group< Flip_map_common, MAX_TESTS_IN_GROUP >;
using Flip_map_object = Flip_map_group::object;

//! Define test group
static Flip_map_group Flip_map( "Base/Flip_map" );

//! Test definitions for group

//! Test tk::flip_pair used to invert a std::pair
template<> template<>
void Flip_map_object::test< 1 >() {
  set_test_name( "flip_pair" );

  std::pair< int, std::string > original{ 1, "blah" };

  ensure( "flipped pair",
    std::pair< std::string, int >{ "blah", 1 } == tk::flip_pair( original ) );
}

//! Test that tk::flip_map does yield a std::multimap sorted by the input
//! std::map::value_type
template<> template<>
void Flip_map_object::test< 2 >() {
  set_test_name( "flip_map" );

  std::map< std::string, tk::real >
    original{ {"two", 2.0}, {"one", 1.0}, {"three", 3.0} };
  std::multimap< tk::real, std::string >
    flipped{ {1.0, "one"}, {2.0, "two"}, {3.0, "three"} };

  auto computed = tk::flip_map( original );

  ensure_not( "flipped map not empty", flipped.empty() );
  ensure( "map flipped", flipped == computed );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
