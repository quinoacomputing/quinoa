// *****************************************************************************
/*!
  \file      src/UnitTest/tests/LoadBalance/TestLinearMap.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for LoadBalance/LinearMap
  \details   Unit tests for LoadBalance/LinearMap
*/
// *****************************************************************************
#ifndef test_LinearMap_h
#define test_LinearMap_h

#include "NoWarning/tut.h"

#include "LinearMap.h"
#include "TestArray.h"

namespace unittest {

extern std::string g_executable;

} // unittest::

namespace tut {

//! All tests in group inherited from this base
struct LinearMap_common {};

//! Test group shortcuts
using LinearMap_group =
  test_group< LinearMap_common, MAX_TESTS_IN_GROUP >;
using LinearMap_object = LinearMap_group::object;

//! Define test group
static LinearMap_group LinearMap( "LoadBalance/LinearMap" );

//! Test definitions for group

//! Test if constructor does not throw on positive number of elements
//! \author J. Bakosi
template<> template<>
void LinearMap_object::test< 1 >() {
  set_test_name( "ctor doesn't throw on positive nelem" );
  tk::CProxy_LinearMap::ckNew( 2 );
}

//! Test use of LinearMap creating an array with nchare <= numpes
//! \author J. Bakosi
template<> template<>
void LinearMap_object::test< 2 >() {
  int nchare = CkNumPes() > 1 ? CkNumPes()/2 : 1;
  set_test_name( "use with nchare (" + std::to_string(nchare) + ") <= numpes ("
                 + std::to_string(CkNumPes()) + ")" );

  // Create linear map chare
  auto map = tk::CProxy_LinearMap::ckNew( nchare );

  // Create array options object for use with linear map chare
  CkArrayOptions opts( nchare );
  opts.setMap( map );

  // Create chare array using linear map
  CProxy_TestArray arrayproxy = CProxy_TestArray::ckNew( opts );
  arrayproxy.doneInserting();

  // If this test fails it will spew errors on the screen...
}

//! Test use of LinearMap creating an array with nchare > numpes
//! \author J. Bakosi
template<> template<>
void LinearMap_object::test< 3 >() {
  int nchare = 2 * CkNumPes();
  set_test_name( "use with nchare (" + std::to_string(nchare) + ") > numpes ("
                 + std::to_string(CkNumPes()) + ")" );

  // Create linear map chare
  auto map = tk::CProxy_LinearMap::ckNew( nchare );

  // Create array options object for use with linear map chare
  CkArrayOptions opts( nchare );
  opts.setMap( map );

  // Create chare array using linear map
  CProxy_TestArray arrayproxy = CProxy_TestArray::ckNew( opts );
  arrayproxy.doneInserting();

  // If this test fails it will spew errors on the screen...
}

} // tut::

#endif // test_LinearMap_h
