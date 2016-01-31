//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestDataLayout.h
  \author    J. Bakosi
  \date      Sun 31 Jan 2016 06:58:42 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit tests for Base/DataLayout.h
  \details   Unit tests for Base/DataLayout.h
*/
//******************************************************************************
#ifndef test_DataLayout_h
#define test_DataLayout_h

#include <tut/tut.hpp>

#include <boost/mpl/at.hpp>

#include "DataLayout.h"

namespace tut {

//! All tests in group inherited from this base
struct DataLayout_common {

  template< typename... Ts >
  struct Test {
    using list = typename tk::make_list< Ts... >::type;
  };

};

//! Test group shortcuts
using DataLayout_group = test_group< DataLayout_common, MAX_TESTS_IN_GROUP >;
using DataLayout_object = DataLayout_group::object;

//! Define test group
DataLayout_group DataLayout( "Base/DataLayout" );

//! Test definitions for group

//! \brief Test that tk::DataLayout' default constructor creates zero
//!   size array
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 1 >() {
  set_test_name( "default ctor yields 0-size array" );

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::size() returns 0",
                 tk::DataLayout< tk::UnkEqComp >().size(), 0 );
  ensure_equals( "<EqCompUnk>::size() returns 0",
                 tk::DataLayout< tk::EqCompUnk >().size(), 0 );
}

//! \brief Test that tk::DataLayout' never leaves its underlying array
//!   pointer uninitialized
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 2 >() {
  set_test_name( "underlying pointer initialized" );

  // Test all template specializations
  ensure_not( "<UnkEqComp>::ptr() not nullptr",
              tk::DataLayout< tk::UnkEqComp >().ptr() == nullptr );
  ensure_not( "<EqCompUnk>::ptr() not nullptr",
              tk::DataLayout< tk::EqCompUnk >().ptr() == nullptr );
}

//! \brief Test that tk::DataLayout' non-default constructor creates
//!   correctly sized arrays
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 3 >() {
  set_test_name( "ctor yields correct-size array" );

  // Test all template specializations and all callable ways
  ensure_equals( "<UnkEqComp>::size() returns 0",
                 tk::DataLayout< tk::UnkEqComp >( 2 ).size(), 0 );
  ensure_equals( "<UnkEqComp>::nunk() returns 2",
                 tk::DataLayout< tk::UnkEqComp >( 2 ).nunk(), 2 );

  ensure_equals( "<UnkEqComp>::size() returns 0",
                 tk::DataLayout< tk::UnkEqComp >( 2, 3 ).size(), 6 );
  ensure_equals( "<UnkEqComp>::nunk() returns 2",
                 tk::DataLayout< tk::UnkEqComp >( 2, 3 ).nunk(), 2 );

  ensure_equals( "<EqCompUnk>::size() returns 0",
                 tk::DataLayout< tk::EqCompUnk >( 2 ).size(), 0 );
  ensure_equals( "<EqCompUnk>::nunk() returns 2",
                 tk::DataLayout< tk::EqCompUnk >( 2 ).nunk(), 2 );

  ensure_equals( "<EqCompUnk>::size() returns 0",
                 tk::DataLayout< tk::EqCompUnk >( 2, 3 ).size(), 6 );
  ensure_equals( "<EqCompUnk>::nunk() returns 2",
                 tk::DataLayout< tk::EqCompUnk >( 2, 3 ).nunk(), 2 );
}

//! \brief Test that tk::DataLayout' operator() returns the correct
//!   value
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 4 >() {
  set_test_name( "operator() returns correct val" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 5 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 5 );

  pp( 1, 2, 2 ) = 0.23;
  pe( 1, 2, 2 ) = 0.32;

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::() value correct", pp(1,2,2), 0.23 );
  ensure_equals( "<EqCompUnk>::() value correct", pe(1,2,2), 0.32 );
}

//! \brief Test that tk::DataLayout' cvar( cptr() ) returns the correct
//!   value
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 5 >() {
  set_test_name( "cvar(cptr()) returns correct val" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 5 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 5 );

  pp( 1, 2, 2 ) = 0.23;
  pe( 1, 2, 2 ) = 0.32;

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::cvar(cptr) value correct",
                 pp.cvar( pp.cptr(2,2), 1 ), 0.23 );
  ensure_equals( "<EqCompUnk>::cvar(cptr) value correct",
                 pe.cvar( pe.cptr(2,2), 1 ), 0.32 );
}

//! \brief Test that tk::DataLayout' cvar( cptr() ) is equivalent to
//!    operator()
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 6 >() {
  set_test_name( "cvar(cptr()) == operator()" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 6 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 6 );

  pp( 1, 2, 3 ) = 0.23;
  pe( 1, 2, 3 ) = 0.32;

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::cvar(cptr) value correct",
                 pp.cvar( pp.cptr(2,3), 1 ), pp(1,2,3) );
  ensure_equals( "<EqCompUnk>::cvar(cptr) value correct",
                 pe.cvar( pe.cptr(2,3), 1 ), pe(1,2,3) );
}

//! \brief Test that tk::DataLayout' major() returns correct string
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 7 >() {
  set_test_name( "major()" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 6 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 6 );

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::major() correct", pp.major(), "unknown-major" );
  ensure_equals( "<EqCompUnk>::major() correct", pe.major(), "equation-major" );
}
} // tut::

#endif // test_DataLayout_h
