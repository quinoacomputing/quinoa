// *****************************************************************************
/*!
  \file      src/UnitTest/tests/LoadBalance/TestUnsMeshMap.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for LoadBalance/UnsMeshMap
  \details   Unit tests for LoadBalance/UnsMeshMap
*/
// *****************************************************************************
#ifndef test_UnsMeshMap_h
#define test_UnsMeshMap_h

#include "NoWarning/tut.h"

#include "UnsMeshMap.h"
#include "TestArray.h"

namespace unittest {

extern std::string g_executable;

} // unittest::

namespace tut {

//! All tests in group inherited from this base
struct UnsMeshMap_common {};

//! Test group shortcuts
using UnsMeshMap_group =
  test_group< UnsMeshMap_common, MAX_TESTS_IN_GROUP >;
using UnsMeshMap_object = UnsMeshMap_group::object;

//! Define test group
static UnsMeshMap_group UnsMeshMap( "LoadBalance/UnsMeshMap" );

//! Test definitions for group

//! Test if constructor behaves on sane input
//! \author J. Bakosi
template<> template<>
void UnsMeshMap_object::test< 1 >() {
  set_test_name( "ctor behaves on sane input" );

  auto nchare = static_cast< std::size_t >( CkNumPes() );
  std::vector< std::vector< std::size_t > > point( nchare );
  std::size_t npoin = nchare;   // each chare will own this many points
  std::size_t i = 0;
  for (auto& p : point)
    for (std::size_t q=0; q<npoin; ++q)
      p.push_back( i++ );

  tk::CProxy_UnsMeshMap::ckNew( nchare*npoin, point );
}

//! \brief Test use of UnsMeshMap creating an array with nchare = numpes when
//!   the points are nicely distributed
//! \author J. Bakosi
template<> template<>
void UnsMeshMap_object::test< 2 >() {
  auto nchare = static_cast< std::size_t >( CkNumPes() );
  set_test_name( "use with nchare (" + std::to_string(nchare) + ") = numpes ("
                 + std::to_string(CkNumPes()) + "), dist" );

  // Construct a vector of vectors storing point ids owned by each chare
  std::vector< std::vector< std::size_t > > point( nchare );
  std::size_t npoin = nchare;   // each chare will own this many points
  std::size_t i = 0;
  for (auto& p : point)
    for (std::size_t q=0; q<npoin; ++q)
      p.push_back( i++ );

  // Create unstructured-mesh-aware map chare
  auto map = tk::CProxy_UnsMeshMap::ckNew( nchare*npoin, point );

  // Create array options object for use with unstructured-mesh-aware map chare
  CkArrayOptions opts( static_cast<int>(nchare) );
  opts.setMap( map );

  // Create chare array using unstructured-mesh-aware map
  CProxy_TestArray arrayproxy = CProxy_TestArray::ckNew( opts );
  arrayproxy.doneInserting();

  // If this test fails it will spew errors on the screen...
}

//! \brief Test use of UnsMeshMap creating an array with nchare > numpes when
//!   the points are nicely distributed
//! \author J. Bakosi
template<> template<>
void UnsMeshMap_object::test< 3 >() {
  auto nchare = static_cast< std::size_t >( 2 * CkNumPes() );
  set_test_name( "use with nchare (" + std::to_string(nchare) + ") > numpes ("
                 + std::to_string(CkNumPes()) + "), dist" );

  // Construct a vector of vectors storing point ids owned by each chare
  std::vector< std::vector< std::size_t > > point( nchare );
  std::size_t npoin = nchare;   // each chare will own this many points
  std::size_t i = 0;
  for (auto& p : point)
    for (std::size_t q=0; q<npoin; ++q)
      p.push_back( i++ );

  // Create unstructured-mesh-aware map chare
  auto map = tk::CProxy_UnsMeshMap::ckNew( nchare*npoin, point );

  // Create array options object for use with unstructured-mesh-aware map chare
  CkArrayOptions opts( static_cast<int>(nchare) );
  opts.setMap( map );

  // Create chare array using unstructured-mesh-aware map
  CProxy_TestArray arrayproxy = CProxy_TestArray::ckNew( opts );
  arrayproxy.doneInserting();

  // If this test fails it will spew errors on the screen...
}

//! \brief Test use of UnsMeshMap creating an array with nchare = numpes when
//!   the points are not distributed at all, but all assigned to PE 0
//! \details Note that the distribution of points this test tests is pretty
//!   pathological since all points are assigned an id 0. Obviously this never
//!   happens in real code, but the algorithm in UnsMeshMap does not care and it
//!   will distribute all points (and array chares) to PE 0, which is the worst
//!   case for what we are testing here, UnsMeshMap::fixPEs(). The test is
//!   successful if the assert at the end of fixPEs() does not fire.
//! \author J. Bakosi
template<> template<>
void UnsMeshMap_object::test< 4 >() {
  auto nchare = static_cast< std::size_t >( CkNumPes() );
  set_test_name( "use with nchare (" + std::to_string(nchare) + ") = numpes ("
                 + std::to_string(CkNumPes()) + "), all 0" );

  // Construct a vector of vectors storing point ids owned by each chare
  std::vector< std::vector< std::size_t > > point( nchare );
  std::size_t npoin = nchare;   // each chare will own this many points
  for (auto& p : point) p.insert( begin(p), npoin, 0 );

  // Create unstructured-mesh-aware map chare
  auto map = tk::CProxy_UnsMeshMap::ckNew( nchare*npoin, point );

  // Create array options object for use with unstructured-mesh-aware map chare
  CkArrayOptions opts( static_cast<int>(nchare) );
  opts.setMap( map );

  // Create chare array using unstructured-mesh-aware map
  CProxy_TestArray arrayproxy = CProxy_TestArray::ckNew( opts );
  arrayproxy.doneInserting();

  // If this test fails it will spew errors on the screen...
}

//! \brief Test use of UnsMeshMap creating an array with nchare > numpes when
//!   the points are not distributed at all, but all assigned to PE 0
//! \details Note that the distribution of points this test tests is pretty
//!   pathological since all points are assigned an id 0. Obviously this never
//!   happens in real code, but the algorithm in UnsMeshMap does not care and it
//!   will distribute all points (and array chares) to PE 0, which is the worst
//!   case for what we are testing here, UnsMeshMap::fixPEs(). The test is
//!   successful if the assert at the end of fixPEs() does not fire.
//! \author J. Bakosi
template<> template<>
void UnsMeshMap_object::test< 5 >() {
  auto nchare = static_cast< std::size_t >( 2 * CkNumPes() );
  set_test_name( "use with nchare (" + std::to_string(nchare) + ") > numpes ("
                 + std::to_string(CkNumPes()) + "), all 0" );

  // Construct a vector of vectors storing point ids owned by each chare
  std::vector< std::vector< std::size_t > > point( nchare );
  std::size_t npoin = nchare;   // each chare will own this many points
  for (auto& p : point) p.insert( begin(p), npoin, 0 );

  // Create unstructured-mesh-aware map chare
  auto map = tk::CProxy_UnsMeshMap::ckNew( nchare*npoin, point );

  // Create array options object for use with unstructured-mesh-aware map chare
  CkArrayOptions opts( static_cast<int>(nchare) );
  opts.setMap( map );

  // Create chare array using unstructured-mesh-aware map
  CProxy_TestArray arrayproxy = CProxy_TestArray::ckNew( opts );
  arrayproxy.doneInserting();

  // If this test fails it will spew errors on the screen...
}

} // tut::

#endif // test_UnsMeshMap_h
