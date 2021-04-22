// *****************************************************************************
/*!
  \file      tests/unit/Base/TestProcessControl.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/ProcessControl.hpp
  \details   Unit tests for Base/ProcessControl.hpp
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Exception.hpp"
#include "ProcessControl.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct ProcessControl_common {};

//! Test group shortcuts
using ProcessControl_group =
  test_group< ProcessControl_common, MAX_TESTS_IN_GROUP >;
using ProcessControl_object = ProcessControl_group::object;

//! Define test group
static ProcessControl_group ProcessControl( "Base/ProcessControl" );

//! Test definitions for group

//! Attempt to call rm with empty argument 
template<> template<>
void ProcessControl_object::test< 1 >() {
  set_test_name( "rm throws on empty argument" );

  try {
    tk::rm( "" );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, tk::rm() still did not fail, test ok
  }
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
