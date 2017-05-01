// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestProcessControl.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/ProcessControl.h
  \details   Unit tests for Base/ProcessControl.h
*/
// *****************************************************************************
#ifndef test_ProcessControl_h
#define test_ProcessControl_h

#include "NoWarning/tut.h"

#include "ProcessControl.h"

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
//! \author J. Bakosi
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

#endif // test_ProcessControl_h
