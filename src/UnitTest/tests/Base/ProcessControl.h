//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/ProcessControl.h
  \author    J. Bakosi
  \date      Tue 24 Mar 2015 02:15:40 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit tests for Base/ProcessControl.h
  \details   Unit tests for Base/ProcessControl.h
*/
//******************************************************************************
#ifndef test_ProcessControl_h
#define test_ProcessControl_h

#include <tut/tut.hpp>
#include <ProcessControl.h>

namespace tut {

//! All tests in group inherited from this base
struct ProcessControl_common {};

//! Test group shortcuts
using ProcessControl_group = test_group< ProcessControl_common >;
using ProcessControl_object = ProcessControl_group::object;

//! Define test group
ProcessControl_group ProcessControl( "Base/ProcessControl" );

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
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, tk::rm() still did not fail, test ok
  }
}

} // tut::

#endif // test_ProcessControl_h
