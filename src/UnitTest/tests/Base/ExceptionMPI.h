//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/ExceptionMPI.h
  \author    J. Bakosi
  \date      Mon 09 Mar 2015 09:13:38 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit tests for Base/ExceptionMPI.h
  \details   Unit tests for Base/ExceptionMPI.h
*/
//******************************************************************************
#ifndef test_ExceptionMPI_h
#define test_ExceptionMPI_h

#include <tut/tut.hpp>
#include <ExceptionMPI.h>

namespace tut {

//! All tests in group inherited from this base
struct ExceptionMPI_common {};

//! Test group shortcuts
using ExceptionMPI_group = test_group< ExceptionMPI_common >;
using ExceptionMPI_object = ExceptionMPI_group::object;

//! Define test group
ExceptionMPI_group ExceptionMPI( "Base/ExceptionMPI" );

//! Test definitions for group

//! Test that AssertMPI macro throws if condition is false on all ranks
template<> template<>
void ExceptionMPI_object::test< 1 >() {
  set_test_name( "AssertMPI macro throws all false" );

  try {
    AssertMPI( 0 == 1, "msg" );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }
}

//! Test that AssertMPI macro does not throw if condition is true on all ranks
template<> template<>
void ExceptionMPI_object::test< 2 >() {
  set_test_name( "AssertMPI macro doesn't throw all true" );

  try {
    AssertMPI( 1 == 1, "msg" );
  }
  catch ( tk::Exception& e ) {
    fail( "should not throw exception" );
  }
}

//! Test that ErrChkMPI macro throws if condition is false on all ranks
template<> template<>
void ExceptionMPI_object::test< 3 >() {
  set_test_name( "ErrChkMPI macro throws all false" );

  try {
    ErrChkMPI( 0 == 1, "msg" );
    fail( "should throw excecption" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown, test ok
  }
}

//! Test that ErrChkMPI macro does not throw if condition is true on all ranks
template<> template<>
void ExceptionMPI_object::test< 4 >() {
  set_test_name( "ErrChkMPI macro doesn't throw all true" );

  try {
    ErrChkMPI( 0 != 1, "msg" );
  }
  catch ( tk::Exception& e ) {
    fail( "should not throw excecption" );
  }
}

//! Test that ErrChkMPI macro throws if condition is false on the 0th rank
template<> template<>
void ExceptionMPI_object::test< 5 >() {
  set_test_name( "ErrChkMPI macro throws 0th false" );

  try {
    int peid;
    MPI_Comm_rank( MPI_COMM_WORLD, &peid );
    ErrChkMPI( peid == 0 ? 0 == 1 : 1 == 1, "msg" );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown, test ok
  }
}

//! Test that ErrChkMPI macro throws if condition is true on only the 0th rank
template<> template<>
void ExceptionMPI_object::test< 6 >() {
  set_test_name( "ErrChkMPI macro throws 0th true only" );

  int numpes;
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );
  if ( numpes > 1 )
    try {
      int peid;
      MPI_Comm_rank( MPI_COMM_WORLD, &peid );
      ErrChkMPI( peid == 0 ? 1 == 1 : 0 == 1, "msg" );
      fail( "should throw exception" );
    }
    catch ( tk::Exception& e ) {
      // exception thrown, test ok
    }
  else
    skip( "in serial, needs multiple PEs" );
}

} // tut::

#endif // test_ExceptionMPI_h
