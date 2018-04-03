// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestData.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/Data.h
  \details   Unit tests for Base/Data.h
*/
// *****************************************************************************
#ifndef test_Data_h
#define test_Data_h

#include <limits>
#include <array>
#include <vector>

#include "NoWarning/tut.h"

#include <boost/mpl/at.hpp>

#include "Data.h"
#include "TUTUtil.h"

namespace tut {

//! All tests in group inherited from this base
struct Data_common {
  const tk::real prec = std::numeric_limits< tk::real >::epsilon();
};

//! Test group shortcuts
using Data_group = test_group< Data_common, MAX_TESTS_IN_GROUP >;
using Data_object = Data_group::object;

//! Define test group
static Data_group Data( "Base/Data" );

//! Test definitions for group

//! Test that tk::Data's constructor creates
//!   correctly sized arrays
template<> template<>
void Data_object::test< 1 >() {
  set_test_name( "correct size" );

  // Test all template specializations and all callable ways
  ensure_equals( "<UnkEqComp>::nunk() returns 2",
                 tk::Data< tk::UnkEqComp >( 2, 3 ).nunk(), 2 );
  ensure_equals( "<UnkEqComp>::nprop() returns 3",
                 tk::Data< tk::UnkEqComp >( 2, 3 ).nprop(), 3 );

  ensure_equals( "<EqCompUnk>::nunk() returns 2",
                 tk::Data< tk::EqCompUnk >( 2, 3 ).nunk(), 2 );
  ensure_equals( "<EqCompUnk>::nprop() returns 3",
                 tk::Data< tk::EqCompUnk >( 2, 3 ).nprop(), 3 );
}

//! Test that tk::Data's operator() returns the correct value
template<> template<>
void Data_object::test< 2 >() {
  set_test_name( "operator() returns correct value" );

  tk::Data< tk::UnkEqComp > pp( 2, 5 );
  tk::Data< tk::EqCompUnk > pe( 2, 5 );

  pp( 1, 2, 2 ) = 0.23;
  pe( 1, 2, 2 ) = 0.32;

  const tk::real& cu = pp( 1, 2, 2 );
  const tk::real& ce = pe( 1, 2, 2 );

  tk::real& u = pp( 1, 2, 2 );
  tk::real& e = pe( 1, 2, 2 );

  // Test all template specializations const-ref access
  ensure_equals( "<UnkEqComp>::() cref incorrect", cu, 0.23 );
  ensure_equals( "<EqCompUnk>::() cref incorrect", ce, 0.32 );
  // Test all template specializations ref access
  ensure_equals( "<UnkEqComp>::() ref incorrect", u, 0.23 );
  ensure_equals( "<EqCompUnk>::() ref incorrect", e, 0.32 );
}

//! Test that tk::Data's operator() throws for out-of-bounds indices
template<> template<>
void Data_object::test< 3 >() {
  set_test_name( "operator() throws for out-of-bounds indices" );

  tk::Data< tk::UnkEqComp > pp( 2, 5 );
  tk::Data< tk::EqCompUnk > pe( 2, 5 );

  pp( 1, 2, 2 ) = 0.23;
  pe( 1, 2, 2 ) = 0.32;

  // Test all template specializations const-ref access

  try {
    pp( 2, 2, 2 );   // unknown out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pp( 1, 3, 2 );   // offset+component out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pp( 1, 2, 3 );   // offset+component out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pe( 2, 2, 2 );   // unknown out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pe( 1, 3, 2 );   // offset+component out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pe( 1, 2, 3 );   // offset+component out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  // Test all template specializations non-const-ref access

  try {
    pp( 2, 2, 2 );  // unknown out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pp( 1, 3, 2 );  // offset+component out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pp( 1, 2, 3 );  // offset+component out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pe( 2, 2, 2 );  // unknown out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pe( 1, 3, 2 );  // offset+component out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    pe( 1, 2, 3 );  // offset+component out of bounds
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }
}

//! Test that tk::Data's var( cptr() ) returns the correct value
template<> template<>
void Data_object::test< 4 >() {
  set_test_name( "var(cptr()) returns correct value" );

  tk::Data< tk::UnkEqComp > pp( 2, 5 );
  tk::Data< tk::EqCompUnk > pe( 2, 5 );

  pp( 1, 2, 2 ) = 0.23;
  pe( 1, 2, 2 ) = 0.32;

  const tk::real& cu = pp.var( pp.cptr(2,2), 1 );
  const tk::real& ce = pe.var( pe.cptr(2,2), 1 );

  tk::real& u = pp.var( pp.cptr(2,2), 1 );
  tk::real& e = pe.var( pe.cptr(2,2), 1 );

  // Test all template specializations const-ref access
  ensure_equals( "<UnkEqComp>::var(cptr) cref incorrect", cu, 0.23 );
  ensure_equals( "<EqCompUnk>::var(cptr) cref incorrect", ce, 0.32 );
  // Test all template specializations non-const-ref access
  ensure_equals( "<UnkEqComp>::var(cptr) ref incorrect", u, 0.23 );
  ensure_equals( "<EqCompUnk>::var(cptr) ref incorrect", e, 0.32 );
}

//! Test that tk::Data's var() throws for out-of-bounds indices
template<> template<>
void Data_object::test< 5 >() {
  set_test_name( "var(cptr()) throws for out-of-bounds indices" );

  tk::Data< tk::UnkEqComp > pp( 2, 5 );
  tk::Data< tk::EqCompUnk > pe( 2, 5 );

  pp( 1, 2, 2 ) = 0.23;
  pe( 1, 2, 2 ) = 0.32;

  // Test all template specializations const-ref access

  try {
    // unknown out of bounds
    pp.var( pp.cptr(2,2), 2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // offset+component out of bounds
    pp.var( pp.cptr(1,3), 2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // offset+component out of bounds
    pp.var( pp.cptr(1,2), 3 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // unknown out of bounds
    pe.var( pe.cptr(2,2), 2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // offset+component out of bounds
    pe.var( pe.cptr(1,3), 2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // offset+component out of bounds
    pe.var( pe.cptr(1,2), 3 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  // Test all template specializations non-const-ref access

  try {
   // unknown out of bounds
   pp.var( pp.cptr(2,2), 2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // offset+component out of bounds
    pp.var( pp.cptr(1,3), 2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // offset+component out of bounds
    pp.var( pp.cptr(1,2), 3 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // unknown out of bounds
    pe.var( pe.cptr(2,2), 2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // offset+component out of bounds
    pe.var( pe.cptr(1,3), 2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  try {
    // offset+component out of bounds
    pe.var( pe.cptr(1,2), 3 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }
}

//! Test that tk::Data's var( cptr() ) is equivalent to operator()
template<> template<>
void Data_object::test< 6 >() {
  set_test_name( "var(cptr()) == operator()" );

  tk::Data< tk::UnkEqComp > pp( 2, 6 );
  tk::Data< tk::EqCompUnk > pe( 2, 6 );

  pp( 1, 2, 3 ) = 0.23;
  pe( 1, 2, 3 ) = 0.32;

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::var(cptr) value incorrect",
                 pp.var( pp.cptr(2,3), 1 ), pp(1,2,3), prec );
  ensure_equals( "<EqCompUnk>::var(cptr) value incorrect",
                 pe.var( pe.cptr(2,3), 1 ), pe(1,2,3), prec );
}

//! Test that tk::Data's layou() returns correct string
template<> template<>
void Data_object::test< 7 >() {
  set_test_name( "layou()" );

  tk::Data< tk::UnkEqComp > pp( 2, 6 );
  tk::Data< tk::EqCompUnk > pe( 2, 6 );

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::layou() correct", pp.layout(), "unknown-major" );
  ensure_equals( "<EqCompUnk>::layou() correct", pe.layout(), "equation-major" );
}

//! Test that tk::Data's extract() returns correct vector of unknowns
template<> template<>
void Data_object::test< 8 >() {
  set_test_name( "extract() vector of unknowns" );

  tk::Data< tk::UnkEqComp > pp( 2, 3 );
  tk::Data< tk::EqCompUnk > pe( 2, 3 );

  pp( 0, 0, 0 ) = 0.1;
  pp( 0, 1, 0 ) = 0.2;
  pp( 0, 2, 0 ) = 0.3;
  pp( 1, 0, 0 ) = 0.4;
  pp( 1, 1, 0 ) = 0.5;
  pp( 1, 2, 0 ) = 0.6;

  pe( 0, 0, 0 ) = 0.1;
  pe( 0, 1, 0 ) = 0.2;
  pe( 0, 2, 0 ) = 0.3;
  pe( 1, 0, 0 ) = 0.4;
  pe( 1, 1, 0 ) = 0.5;
  pe( 1, 2, 0 ) = 0.6;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::extract() vector of unknowns at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.4 }, pp.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::extract() vector of unknowns at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.5 }, pp.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::extract() vector of unknowns at 2,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.6 }, pp.extract( 2, 0 ) );
  veceq( "<UnkEqComp>::extract() vector of unknowns at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.5 }, pp.extract( 0, 1 ) );
  veceq( "<UnkEqComp>::extract() vector of unknowns at 1,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.6 }, pp.extract( 1, 1 ) );
  veceq( "<UnkEqComp>::extract() vector of unknowns at 0,2 incorrect",
         std::vector< tk::real >{ 0.3, 0.6 }, pp.extract( 0, 2 ) );

  veceq( "<EqCompUnk>::extract() vector of unknowns at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.4 }, pe.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::extract() vector of unknowns at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.5 }, pe.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::extract() vector of unknowns at 2,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.6 }, pe.extract( 2, 0 ) );
  veceq( "<EqCompUnk>::extract() vector of unknowns at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.5 }, pe.extract( 0, 1 ) );
  veceq( "<EqCompUnk>::extract() vector of unknowns at 1,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.6 }, pe.extract( 1, 1 ) );
  veceq( "<EqCompUnk>::extract() vector of unknowns at 0,2 incorrect",
         std::vector< tk::real >{ 0.3, 0.6 }, pe.extract( 0, 2 ) );
}

//! Test that tk::Data's extract() returns correct vector components
template<> template<>
void Data_object::test< 9 >() {
  set_test_name( "extract() vector of components" );

  tk::Data< tk::UnkEqComp > pp( 2, 3 );
  tk::Data< tk::EqCompUnk > pe( 2, 3 );

  pp( 0, 0, 0 ) = 0.1;
  pp( 0, 1, 0 ) = 0.2;
  pp( 0, 2, 0 ) = 0.3;
  pp( 1, 0, 0 ) = 0.4;
  pp( 1, 1, 0 ) = 0.5;
  pp( 1, 2, 0 ) = 0.6;

  pe( 0, 0, 0 ) = 0.3;
  pe( 0, 1, 0 ) = 0.2;
  pe( 0, 2, 0 ) = 0.1;
  pe( 1, 0, 0 ) = 0.6;
  pe( 1, 1, 0 ) = 0.5;
  pe( 1, 2, 0 ) = 0.4;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::extract() vector of components at 0 incorrect",
         std::vector< tk::real >{ 0.1, 0.2, 0.3 }, pp.extract( 0 ) );
  veceq( "<UnkEqComp>::extract() vector of components at 1 incorrect",
         std::vector< tk::real >{ 0.4, 0.5, 0.6 }, pp.extract( 1 ) );

  veceq( "<EqCompUnk>::extract() vector of components at 0 incorrect",
         std::vector< tk::real >{ 0.3, 0.2, 0.1 }, pe.extract( 0 ) );
  veceq( "<EqCompUnk>::extract() vector of components at 1 incorrect",
         std::vector< tk::real >{ 0.6, 0.5, 0.4 }, pe.extract( 1 ) );
}

//! Test that tk::Data's operator[] returns correct vector components
template<> template<>
void Data_object::test< 10 >() {
  set_test_name( "operator[] to access vector of components" );

  tk::Data< tk::UnkEqComp > pp( 2, 3 );
  tk::Data< tk::EqCompUnk > pe( 2, 3 );

  pp( 0, 0, 0 ) = 0.1;
  pp( 0, 1, 0 ) = 0.2;
  pp( 0, 2, 0 ) = 0.3;
  pp( 1, 0, 0 ) = 0.4;
  pp( 1, 1, 0 ) = 0.5;
  pp( 1, 2, 0 ) = 0.6;

  pe( 0, 0, 0 ) = 0.3;
  pe( 0, 1, 0 ) = 0.2;
  pe( 0, 2, 0 ) = 0.1;
  pe( 1, 0, 0 ) = 0.6;
  pe( 1, 1, 0 ) = 0.5;
  pe( 1, 2, 0 ) = 0.4;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator[] returning vector of components at 0 incorrect",
         std::vector< tk::real >{ 0.1, 0.2, 0.3 }, pp[0] );
  veceq( "<UnkEqComp>::operator[] returning vector of components at 1 incorrect",
         std::vector< tk::real >{ 0.4, 0.5, 0.6 }, pp[1] );

  veceq( "<EqCompUnk>::operator[] returning vector of components at 0 incorrect",
         std::vector< tk::real >{ 0.3, 0.2, 0.1 }, pe[0] );
  veceq( "<EqCompUnk>::operator[] returning vector of components at 1 incorrect",
         std::vector< tk::real >{ 0.6, 0.5, 0.4 }, pe[1] );
}

//! Test that tk::Data's extract() returns correct array of four reals
template<> template<>
void Data_object::test< 11 >() {
  set_test_name( "extract() array of four reals for A,B,C,D" );

  tk::Data< tk::UnkEqComp > pp( 8, 2 );
  tk::Data< tk::EqCompUnk > pe( 8, 2 );

  pp( 0, 0, 0 ) = 0.1;
  pp( 1, 0, 0 ) = 0.2;
  pp( 2, 0, 0 ) = 0.3;
  pp( 3, 0, 0 ) = 0.4;
  pp( 4, 0, 0 ) = 0.5;
  pp( 5, 0, 0 ) = 0.6;
  pp( 6, 0, 0 ) = 0.7;
  pp( 7, 0, 0 ) = 0.8;

  pp( 0, 1, 0 ) = 1.1;
  pp( 1, 1, 0 ) = 1.2;
  pp( 2, 1, 0 ) = 1.3;
  pp( 3, 1, 0 ) = 1.4;
  pp( 4, 1, 0 ) = 1.5;
  pp( 5, 1, 0 ) = 1.6;
  pp( 6, 1, 0 ) = 1.7;
  pp( 7, 1, 0 ) = 1.8;

  pe( 0, 0, 0 ) = 0.1;
  pe( 1, 0, 0 ) = 0.2;
  pe( 2, 0, 0 ) = 0.3;
  pe( 3, 0, 0 ) = 0.4;
  pe( 4, 0, 0 ) = 0.5;
  pe( 5, 0, 0 ) = 0.6;
  pe( 6, 0, 0 ) = 0.7;
  pe( 7, 0, 0 ) = 0.8;

  pe( 0, 1, 0 ) = 1.1;
  pe( 1, 1, 0 ) = 1.2;
  pe( 2, 1, 0 ) = 1.3;
  pe( 3, 1, 0 ) = 1.4;
  pe( 4, 1, 0 ) = 1.5;
  pe( 5, 1, 0 ) = 1.6;
  pe( 6, 1, 0 ) = 1.7;
  pe( 7, 1, 0 ) = 1.8;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::extract() array of four reals at 0,0:3,2,1,0 incorrect",
         std::array< tk::real, 4 >{{ 0.4, 0.3, 0.2, 0.1 }},
           pp.extract( 0, 0, 3, 2, 1, 0 ) );
  veceq( "<UnkEqComp>::extract() array of four reals at 0,0:1,3,5,7 incorrect",
         std::array< tk::real, 4 >{{ 0.2, 0.4, 0.6, 0.8 }},
           pp.extract( 0, 0, 1, 3, 5, 7 ) );
  veceq( "<UnkEqComp>::extract() array of four reals at 0,1:3,2,1,0 incorrect",
         std::array< tk::real, 4 >{{ 1.4, 1.3, 1.2, 1.1 }},
           pp.extract( 0, 1, 3, 2, 1, 0 ) );
  veceq( "<UnkEqComp>::extract() array of four reals at 0,1:1,3,5,7 incorrect",
         std::array< tk::real, 4 >{{ 1.2, 1.4, 1.6, 1.8 }},
           pp.extract( 0, 1, 1, 3, 5, 7 ) );

  veceq( "<EqCompUnk>::extract() array of four reals at 0,0:3,2,1,0 incorrect",
         std::array< tk::real, 4 >{{ 0.4, 0.3, 0.2, 0.1 }},
           pe.extract( 0, 0, 3, 2, 1, 0 ) );
  veceq( "<EqCompUnk>::extract() array of four reals at 0,0:1,3,5,7 incorrect",
         std::array< tk::real, 4 >{{ 0.2, 0.4, 0.6, 0.8 }},
           pe.extract( 0, 0, 1, 3, 5, 7 ) );
  veceq( "<EqCompUnk>::extract() array of four reals at 0,1:3,2,1,0 incorrect",
         std::array< tk::real, 4 >{{ 1.4, 1.3, 1.2, 1.1 }},
           pe.extract( 0, 1, 3, 2, 1, 0 ) );
  veceq( "<kEqCompUnk>::extract() array of four reals 0,1:1,3,5,7 incorrect",
         std::array< tk::real, 4 >{{ 1.2, 1.4, 1.6, 1.8 }},
           pe.extract( 0, 1, 1, 3, 5, 7 ) );
}

//! Test that tk::Data's extract() returns correct array of four reals
template<> template<>
void Data_object::test< 12 >() {
  set_test_name( "extract() array of four reals for N[4]" );

  tk::Data< tk::UnkEqComp > pp( 8, 2 );
  tk::Data< tk::EqCompUnk > pe( 8, 2 );

  pp( 0, 0, 0 ) = 0.1;
  pp( 1, 0, 0 ) = 0.2;
  pp( 2, 0, 0 ) = 0.3;
  pp( 3, 0, 0 ) = 0.4;
  pp( 4, 0, 0 ) = 0.5;
  pp( 5, 0, 0 ) = 0.6;
  pp( 6, 0, 0 ) = 0.7;
  pp( 7, 0, 0 ) = 0.8;

  pp( 0, 1, 0 ) = 1.1;
  pp( 1, 1, 0 ) = 1.2;
  pp( 2, 1, 0 ) = 1.3;
  pp( 3, 1, 0 ) = 1.4;
  pp( 4, 1, 0 ) = 1.5;
  pp( 5, 1, 0 ) = 1.6;
  pp( 6, 1, 0 ) = 1.7;
  pp( 7, 1, 0 ) = 1.8;

  pe( 0, 0, 0 ) = 0.1;
  pe( 1, 0, 0 ) = 0.2;
  pe( 2, 0, 0 ) = 0.3;
  pe( 3, 0, 0 ) = 0.4;
  pe( 4, 0, 0 ) = 0.5;
  pe( 5, 0, 0 ) = 0.6;
  pe( 6, 0, 0 ) = 0.7;
  pe( 7, 0, 0 ) = 0.8;

  pe( 0, 1, 0 ) = 1.1;
  pe( 1, 1, 0 ) = 1.2;
  pe( 2, 1, 0 ) = 1.3;
  pe( 3, 1, 0 ) = 1.4;
  pe( 4, 1, 0 ) = 1.5;
  pe( 5, 1, 0 ) = 1.6;
  pe( 6, 1, 0 ) = 1.7;
  pe( 7, 1, 0 ) = 1.8;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::extract() array of four reals at 0,0:3,2,1,0 incorrect",
          std::array< tk::real, 4 >{{ 0.4, 0.3, 0.2, 0.1 }},
            pp.extract( 0, 0, {{3,2,1,0}} ) );
  veceq( "<UnkEqComp>::extract() array of four reals at 0,0:1,3,5,7 incorrect",
          std::array< tk::real, 4 >{{ 0.2, 0.4, 0.6, 0.8 }},
            pp.extract( 0, 0, {{1,3,5,7}} ) );
  veceq( "<UnkEqComp>::extract() array of four reals at 0,1:3,2,1,0 incorrect",
          std::array< tk::real, 4 >{{ 1.4, 1.3, 1.2, 1.1 }},
            pp.extract( 0, 1, {{3,2,1,0}} ) );
  veceq( "<UnkEqComp>::extract() array of four reals at 0,1:1,3,5,7 incorrect",
          std::array< tk::real, 4 >{{ 1.2, 1.4, 1.6, 1.8 }},
            pp.extract( 0, 1, {{1,3,5,7}} ) );

  veceq( "<EqCompUnk>::extract() array of four reals at 0,0:3,2,1,0 incorrect",
          std::array< tk::real, 4 >{{ 0.4, 0.3, 0.2, 0.1 }},
            pe.extract( 0, 0, {{3,2,1,0}} ) );
  veceq( "<EqCompUnk>::extract() array of four reals at 0,0:1,3,5,7 incorrect",
          std::array< tk::real, 4 >{{ 0.2, 0.4, 0.6, 0.8 }},
            pe.extract( 0, 0, {{1,3,5,7}} ) );
  veceq( "<EqCompUnk>::extract() array of four reals at 0,1:3,2,1,0 incorrect",
          std::array< tk::real, 4 >{{ 1.4, 1.3, 1.2, 1.1 }},
            pe.extract( 0, 1, {{3,2,1,0}} ) );
  veceq( "<kEqCompUnk>::extract() array of four reals 0,1:1,3,5,7 incorrect",
          std::array< tk::real, 4 >{{ 1.2, 1.4, 1.6, 1.8 }},
            pe.extract( 0, 1, {{1,3,5,7}} ) );
}

//! Test that tk::Data's fill() correctly fills complete data array
template<> template<>
void Data_object::test< 13 >() {
  set_test_name( "fill() all with the same value" );

  tk::Data< tk::UnkEqComp > pp( 3, 2 );
  tk::Data< tk::EqCompUnk > pe( 3, 2 );

  pp.fill( 0.0 );
  pe.fill( 0.1 );

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::fill() all with the same value at 0,0 incorrect",
          std::vector< tk::real >{ 0.0, 0.0, 0.0 }, pp.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::fill() all with the same value at 1,0 incorrect",
          std::vector< tk::real >{ 0.0, 0.0, 0.0 }, pp.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::fill() all with the same value at 0,1 incorrect",
          std::vector< tk::real >{ 0.0, 0.0, 0.0 }, pp.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::fill() all with the same value at 0,0 incorrect",
          std::vector< tk::real >{ 0.1, 0.1, 0.1 }, pe.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::fill() all with the same value at 1,0 incorrect",
          std::vector< tk::real >{ 0.1, 0.1, 0.1 }, pe.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::fill() all with the same value at 0,1 incorrect",
          std::vector< tk::real >{ 0.1, 0.1, 0.1 }, pe.extract( 0, 1 ) );
}

//! Test that tk::Data's fill() correctly fills vector of unknowns
template<> template<>
void Data_object::test< 14 >() {
  set_test_name( "fill() vector of unknowns with the same value" );

  tk::Data< tk::UnkEqComp > pp( 3, 2 );
  tk::Data< tk::EqCompUnk > pe( 3, 2 );

  pp.fill( 0, 0, 1.5 );
  pp.fill( 1, 0, 2.5 );

  pe.fill( 0, 0, 0.5 );
  pe.fill( 1, 0, -0.5 );

  using unittest::veceq;

  // Test all template specializations
  veceq(
    "<UnkEqComp>::fill() vector of unknowns with the same value at 0,0 incorrect",
    std::vector< tk::real >{ 1.5, 1.5, 1.5 }, pp.extract( 0, 0 ) );
  veceq(
    "<UnkEqComp>::fill() vector of unknowns with the same value at 1,0 incorrect",
    std::vector< tk::real >{ 2.5, 2.5, 2.5 }, pp.extract( 1, 0 ) );
  veceq(
    "<UnkEqComp>::fill() vector of unknowns with the same value at 0,1 incorrect",
    std::vector< tk::real >{ 2.5, 2.5, 2.5 }, pp.extract( 0, 1 ) );

  veceq(
    "<EqCompUnk>::fill() vector of unknowns with the same value at 0,0 incorrect",
    std::vector< tk::real >{ 0.5, 0.5, 0.5 }, pe.extract( 0, 0 ) );
  veceq(
    "<EqCompUnk>::fill() vector of unknowns with the same value at 1,0 incorrect",
    std::vector< tk::real >{ -0.5, -0.5, -0.5 }, pe.extract( 1, 0 ) );
  veceq(
    "<EqCompUnk>::fill() vector of unknowns with the same value at 0,1 incorrect",
    std::vector< tk::real >{ -0.5, -0.5, -0.5 }, pe.extract( 0, 1 ) );
}

//! \brief Test that tk::Data's memory layout, i.e., stores data with
//!   correct strides
template<> template<>
void Data_object::test< 15 >() {
  set_test_name( "strides" );

  tk::Data< tk::UnkEqComp > pp( 2, 5 );
  tk::Data< tk::EqCompUnk > pe( 2, 5 );

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::component stride incorrect",
                 pp.cptr(2,2) - pp.cptr(1,2), 1 );
  ensure_equals( "<UnkEqComp>::offset stride incorrect",
                 pp.cptr(2,2) - pp.cptr(2,1), 1 );

  ensure_equals( "<EqCompUnk>::component stride incorrect",
                 pe.cptr(2,2) - pe.cptr(1,2), 2 );
  ensure_equals( "<EqCompUnk>::offset stride incorrect",
                 pe.cptr(2,2) - pe.cptr(2,1), 2 );
}

//! Test tk::Data's copy constructor
template<> template<>
void Data_object::test< 16 >() {
  set_test_name( "copy constructor" );

  tk::Data< tk::UnkEqComp > p( 3, 2 );
  tk::Data< tk::EqCompUnk > e( 3, 2 );

  p.fill( 0.1 );
  e.fill( 0.2 );

  std::vector< tk::Data< tk::UnkEqComp > > v;
  v.push_back( p );
  std::vector< tk::Data< tk::EqCompUnk > > w;
  w.push_back( e );

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::ctor() at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, v[0].extract( 0, 0 ) );
  veceq( "<UnkEqComp>::ctor() at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, v[0].extract( 1, 0 ) );
  veceq( "<UnkEqComp>::ctor() at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, v[0].extract( 0, 1 ) );

  veceq( "<EqCompUnk>::ctor() at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, w[0].extract( 0, 0 ) );
  veceq( "<EqCompUnk>::ctor() at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, w[0].extract( 1, 0 ) );
  veceq( "<EqCompUnk>::ctor() at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, w[0].extract( 0, 1 ) );
}

//! Test tk::Data's copy assignment
template<> template<>
void Data_object::test< 17 >() {
  set_test_name( "copy assignment" );

  tk::Data< tk::UnkEqComp > p( 3, 2 );
  tk::Data< tk::EqCompUnk > e( 3, 2 );

  p.fill( 0.1 );
  e.fill( 0.2 );

  tk::Data< tk::UnkEqComp > p1( 3, 2 );
  p1 = p;

  tk::Data< tk::EqCompUnk > e1( 3, 2 );
  e1 = e;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::cass() at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::cass() at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::cass() at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::cass() at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::cass() at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::cass() at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 0, 1 ) );
}

//! Test tk::Data's move constructor
template<> template<>
void Data_object::test< 18 >() {
  set_test_name( "move constructor" );

  tk::Data< tk::UnkEqComp > p( 3, 2 );
  tk::Data< tk::EqCompUnk > e( 3, 2 );

  p.fill( 0.1 );
  e.fill( 0.2 );

  std::vector< tk::Data< tk::UnkEqComp > > v;
  v.emplace_back( p );
  std::vector< tk::Data< tk::EqCompUnk > > w;
  w.emplace_back( e );

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::mctor() at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, v[0].extract( 0, 0 ) );
  veceq( "<UnkEqComp>::mctor() at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, v[0].extract( 1, 0 ) );
  veceq( "<UnkEqComp>::mctor() at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, v[0].extract( 0, 1 ) );

  veceq( "<EqCompUnk>::mctor() at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, w[0].extract( 0, 0 ) );
  veceq( "<EqCompUnk>::mctor() at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, w[0].extract( 1, 0 ) );
  veceq( "<EqCompUnk>::mctor() at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, w[0].extract( 0, 1 ) );
}

//! Test tk::Data's move assignment
template<> template<>
void Data_object::test< 19 >() {
  set_test_name( "move assignment" );

  tk::Data< tk::UnkEqComp > p( 3, 2 );
  tk::Data< tk::EqCompUnk > e( 3, 2 );

  p.fill( 0.1 );
  e.fill( 0.2 );

  auto p1 = std::move( p );
  auto e1 = std::move( e );

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::mass() at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::mass() at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::mass() at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::mass() at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::mass() at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::mass() at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 0, 1 ) );
}

//! Test tk::Data's operator-=
template<> template<>
void Data_object::test< 20 >() {
  set_test_name( "operator-=" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  p1 -= p2;
  e1 -= e2;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator-=() at 0,0 incorrect",
         std::vector< tk::real >{ -0.2, -0.2, -0.2 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator-=() at 1,0 incorrect",
         std::vector< tk::real >{ -0.2, -0.2, -0.2 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator-=() at 0,1 incorrect",
         std::vector< tk::real >{ -0.2, -0.2, -0.2 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator-=() at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator-=() at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator-=() at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e1.extract( 0, 1 ) );

  // Test the rhs of the subtract stay the same for all template specializations
  veceq( "<UnkEqComp>::operator-=() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator-=() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator-=() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator-=() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator-=() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator-=() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 1 ) );
}

//! Test tk::Data's operator-
template<> template<>
void Data_object::test< 21 >() {
  set_test_name( "operator-" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  auto p = p1 - p2;
  auto e = e1 - e2;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator-() at 0,0 incorrect",
         std::vector< tk::real >{ -0.2, -0.2, -0.2 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator-() at 1,0 incorrect",
         std::vector< tk::real >{ -0.2, -0.2, -0.2 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator-() at 0,1 incorrect",
         std::vector< tk::real >{ -0.2, -0.2, -0.2 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator-() at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator-() at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator-() at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, e.extract( 0, 1 ) );

  // Test the lhs of the subtract stay the same for all template specializations
  veceq( "<UnkEqComp>::operator-() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator-() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator-() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator-() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator-() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator-() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 1 ) );

  // Test the rhs of the subtract stay the same for all template specializations
  veceq( "<UnkEqComp>::operator-() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator-() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator-() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator-() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator-() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator-() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 1 ) );
}

//! Test tk::Data's operator+=
template<> template<>
void Data_object::test< 22 >() {
  set_test_name( "operator+=" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  p1 += p2;
  e1 += e2;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator+=() at 0,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator+=() at 1,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator+=() at 0,1 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator+=() at 0,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator+=() at 1,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator+=() at 0,1 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e1.extract( 0, 1 ) );

  // Test the rhs of the addition stay the same for all template specializations
  veceq( "<UnkEqComp>::operator+=() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator+=() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator+=() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator+=() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator+=() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator+=() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 1 ) );
}

//! Test tk::Data's operator+
template<> template<>
void Data_object::test< 23 >() {
  set_test_name( "operator+" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  auto p = p1 + p2;
  auto e = e1 + e2;

  using unittest::veceq;

  // Test the result of the addition for all template specializations
  veceq( "<UnkEqComp>::operator+() res at 0,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator+() res at 1,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator+() res at 0,1 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator+() res at 0,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator+() res at 1,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator+() res at 0,1 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e.extract( 0, 1 ) );

  // Test the lhs of the addition stay the same for all template specializations
  veceq( "<UnkEqComp>::operator+() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator+() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator+() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator+() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator+() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator+() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 1 ) );

  // Test the rhs of the addition stay the same for all template specializations
  veceq( "<UnkEqComp>::operator+() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator+() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator+() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator+() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator+() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator+() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 1 ) );
}

//! Test tk::Data's operator*= by Data as rhs
template<> template<>
void Data_object::test< 24 >() {
  set_test_name( "operator*= by Data right" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  p1 *= p2;
  e1 *= e2;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator*=() at 0,0 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*=() at 1,0 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*=() at 0,1 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*=() at 0,0 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*=() at 1,0 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*=() at 0,1 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, e1.extract( 0, 1 ) );

  // Test the rhs of the multiply stay the same for all template specializations
  veceq( "<UnkEqComp>::operator*=() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*=() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*=() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*=() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*=() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*=() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 1 ) );
}

//! Test tk::Data's operator* by Data as rhs
template<> template<>
void Data_object::test< 25 >() {
  set_test_name( "operator* by Data right" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  auto p = p1 * p2;
  auto e = e1 * e2;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator*() at 0,0 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*() at 1,0 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*() at 0,1 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*() at 0,0 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*() at 1,0 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*() at 0,1 incorrect",
         std::vector< tk::real >{ 0.03, 0.03, 0.03 }, e.extract( 0, 1 ) );

  // Test the lhs of the multiply stay the same for all template specializations
  veceq( "<UnkEqComp>::operator*() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 1 ) );

  // Test the rhs of the multiply stay the same for all template specializations
  veceq( "<UnkEqComp>::operator*() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p2.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e2.extract( 0, 1 ) );
}

//! Test tk::Data's operator*= by tk::real as rhs
template<> template<>
void Data_object::test< 26 >() {
  set_test_name( "operator*= by tk::real right" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.4 );
  e1.fill( 0.3 );

  p1 *= 0.5;
  e1 *= 4;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator*=() at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*=() at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*=() at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*=() at 0,0 incorrect",
         std::vector< tk::real >{ 1.2, 1.2, 1.2 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*=() at 1,0 incorrect",
         std::vector< tk::real >{ 1.2, 1.2, 1.2 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*=() at 0,1 incorrect",
         std::vector< tk::real >{ 1.2, 1.2, 1.2 }, e1.extract( 0, 1 ) );
}

//! Test tk::Data's operator* by tk::real as rhs
template<> template<>
void Data_object::test< 27 >() {
  set_test_name( "operator* by tk::real right" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 );

  p1.fill( 0.1 );
  e1.fill( 0.3 );

  auto p = p1 * 0.2;
  auto e = e1 * 0.3;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator*() at 0,0 incorrect",
         std::vector< tk::real >{ 0.02, 0.02, 0.02 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*() at 1,0 incorrect",
         std::vector< tk::real >{ 0.02, 0.02, 0.02 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*() at 0,1 incorrect",
         std::vector< tk::real >{ 0.02, 0.02, 0.02 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*() at 0,0 incorrect",
         std::vector< tk::real >{ 0.09, 0.09, 0.09 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*() at 1,0 incorrect",
         std::vector< tk::real >{ 0.09, 0.09, 0.09 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*() at 0,1 incorrect",
         std::vector< tk::real >{ 0.09, 0.09, 0.09 }, e.extract( 0, 1 ) );

  // Test the lhs of the multiply stay the same for all template specializations
  veceq( "<UnkEqComp>::operator*() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 1 ) );
}


//! Test tk::Data's operator/=
template<> template<>
void Data_object::test< 28 >() {
  set_test_name( "operator/=" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.2 );
  e1.fill( 0.3 );       e2.fill( 0.6 );

  p1 /= p2;
  e1 /= e2;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator/=() at 0,0 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator/=() at 1,0 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator/=() at 0,1 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator/=() at 0,0 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator/=() at 1,0 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator/=() at 0,1 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, e1.extract( 0, 1 ) );

  // Test the rhs of the division stay the same for all template specializations
  veceq( "<UnkEqComp>::operator/=() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p2.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator/=() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p2.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator/=() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p2.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator/=() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.6, 0.6, 0.6 }, e2.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator/=() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.6, 0.6, 0.6 }, e2.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator/=() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.6, 0.6, 0.6 }, e2.extract( 0, 1 ) );
}

//! Test tk::Data's operator/
template<> template<>
void Data_object::test< 29 >() {
  set_test_name( "operator/" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.2 );
  e1.fill( 0.3 );       e2.fill( 0.6 );

  auto p = p1 / p2;
  auto e = e1 / e2;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator/() at 0,0 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator/() at 1,0 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator/() at 0,1 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator/() at 0,0 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator/() at 1,0 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator/() at 0,1 incorrect",
         std::vector< tk::real >{ 0.5, 0.5, 0.5 }, e.extract( 0, 1 ) );

  // Test the lhs of the division stay the same for all template specializations
  veceq( "<UnkEqComp>::operator/() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator/() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator/() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator/() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator/() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator/() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 1 ) );

  // Test the rhs of the division stay the same for all template specializations
  veceq( "<UnkEqComp>::operator/() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p2.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator/() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p2.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator/() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.2, 0.2, 0.2 }, p2.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator/() rhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.6, 0.6, 0.6 }, e2.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator/() rhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.6, 0.6, 0.6 }, e2.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator/() rhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.6, 0.6, 0.6 }, e2.extract( 0, 1 ) );
}

//! Test tk::Data's operator* by tk::real as lhs
template<> template<>
void Data_object::test< 31 >() {
  set_test_name( "operator* by tk::real left" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 );

  p1.fill( 0.1 );
  e1.fill( 0.3 );

  auto p = 0.2 * p1;
  auto e = 0.3 * e1;

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator*() at 0,0 incorrect",
         std::vector< tk::real >{ 0.02, 0.02, 0.02 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*() at 1,0 incorrect",
         std::vector< tk::real >{ 0.02, 0.02, 0.02 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*() at 0,1 incorrect",
         std::vector< tk::real >{ 0.02, 0.02, 0.02 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*() at 0,0 incorrect",
         std::vector< tk::real >{ 0.09, 0.09, 0.09 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*() at 1,0 incorrect",
         std::vector< tk::real >{ 0.09, 0.09, 0.09 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*() at 0,1 incorrect",
         std::vector< tk::real >{ 0.09, 0.09, 0.09 }, e.extract( 0, 1 ) );

  // Test the lhs of the multiply stay the same for all template specializations
  veceq( "<UnkEqComp>::operator*() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator*() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator*() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p1.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator*() lhs at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator*() lhs at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator*() lhs at 0,1 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e1.extract( 0, 1 ) );
}

//! Test tk::Data's operator min
template<> template<>
void Data_object::test< 32 >() {
  set_test_name( "operator min" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  auto p = tk::min( p1, p2 );
  auto e = tk::min( e1, e2 );

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator min() at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator min() at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator min() at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator min() at 0,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator min() at 1,0 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator min() at 0,1 incorrect",
         std::vector< tk::real >{ 0.1, 0.1, 0.1 }, e.extract( 0, 1 ) );
}

//! Test tk::Data's operator max
template<> template<>
void Data_object::test< 33 >() {
  set_test_name( "operator max" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  auto p = tk::max( p1, p2 );
  auto e = tk::max( e1, e2 );

  using unittest::veceq;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator max() at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator max() at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator max() at 0.3 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator max() at 0,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator max() at 1,0 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator max() at 0.3 incorrect",
         std::vector< tk::real >{ 0.3, 0.3, 0.3 }, e.extract( 0, 1 ) );
}

//! Test tk::Data's operator==
template<> template<>
void Data_object::test< 34 >() {
  set_test_name( "operator==" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.1 );
  e1.fill( 0.3 );       e2.fill( 0.2 );

  ensure_equals( "<UnkEqComp>::operator ==() incorrect", p1==p2, true );
  ensure_equals( "<EqCompUnk>::operator ==() incorrect", e1==e2, false );

  p1.fill( 0.1 );       p2.fill( 0.2 );
  e1.fill( 0.3 );       e2.fill( 0.3 );

  ensure_equals( "<UnkEqComp>::operator ==() incorrect", p1==p2, false );
  ensure_equals( "<EqCompUnk>::operator ==() incorrect", e1==e2, true );
}

//! Test tk::Data's operator!=
template<> template<>
void Data_object::test< 35 >() {
  set_test_name( "operator!=" );

  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.1 );
  e1.fill( 0.3 );       e2.fill( 0.2 );

  ensure_equals( "<UnkEqComp>::operator !=() incorrect", p1!=p2, false );
  ensure_equals( "<EqCompUnk>::operator !=() incorrect", e1!=e2, true );

  p1.fill( 0.1 );       p2.fill( 0.2 );
  e1.fill( 0.3 );       e2.fill( 0.3 );

  ensure_equals( "<UnkEqComp>::operator !=() incorrect", p1!=p2, true );
  ensure_equals( "<EqCompUnk>::operator !=() incorrect", e1!=e2, false );
}

//! Test maxdiff between two tk::Data objects
template<> template<>
void Data_object::test< 36 >() {
  set_test_name( "maxdiff" );

  // Test equal objects with UnkEqComp data layout
  tk::Data< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  p1.fill( 0.1 );
  p2.fill( 0.1 );
  auto m = tk::maxdiff(p1,p2);
  ensure_equals( "<UnkEqComp>::maxdiff eq pos incorrect", m.first, 0, prec );
  ensure_equals( "<UnkEqComp>::maxdiff eq dif incorrect", m.second, 0.0, prec );

  // Test equal objects with EqCompUnk data layout
  tk::Data< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );
  e1.fill( 0.3 );
  e2.fill( 0.3 );
  m = tk::maxdiff(e1,e2);
  ensure_equals( "<EqCompUnk>::maxdiff eq pos incorrect", m.first, 0, prec );
  ensure_equals( "<EqCompUnk>::maxdiff eq dif incorrect", m.second, 0.0, prec );

  // Test unequal objects with UnkEqComp data layout
  tk::Data< tk::UnkEqComp > q1( 3, 2 ), q2( 3, 2 );
  q1.fill( 0.1 );
  q2.fill( 0.0 );  q2(2,0,0) = 1.0;
  m = tk::maxdiff(q1,q2);
  ensure_equals( "<UnkEqComp>::maxdiff neq pos incorrect", m.first, 4, prec );
  ensure_equals( "<UnkEqComp>::maxdiff neq dif incorrect", m.second, 0.9,prec );

  // Test unequal objects with EqCompUnk data layout
  tk::Data< tk::UnkEqComp > f1( 3, 2 ), f2( 3, 2 );
  f1.fill( 0.1 );
  f2.fill( 0.0 );  f2(2,1,0) = -1.0;
  m = tk::maxdiff(f1,f2);
  ensure_equals( "<EqCompUnk>::maxdiff neq pos incorrect", m.first, 5, prec );
  ensure_equals( "<EqCompUnk>::maxdiff neq dif incorrect", m.second, 1.1,prec );
}

//! Test tk::Data::push_back()
template<> template<>
void Data_object::test< 37 >() {
  set_test_name( "push_back" );

  // Test with UnkEqComp data layout
  tk::Data< tk::UnkEqComp > p( 3, 2 );
  p(0,0,0) = 1.0;  p(0,1,0) = 2.0;
  p(1,0,0) = 3.0;  p(1,1,0) = 4.0;
  p(2,0,0) = 5.0;  p(2,1,0) = 6.0;

  p.push_back( {0.2, 0.3} );

  ensure_equals( "nunk after <UnkEqComp>::push_back() incorrect", p.nunk(), 4 );
  ensure_equals( "nprop after <UnkEqComp>::push_back() incorrect", p.nprop(), 2 );

  using unittest::veceq;

  veceq( "<UnkEqComp>::push_back() at 0 incorrect",
         std::vector< tk::real >{ 1.0, 2.0 }, p[0] );
  veceq( "<UnkEqComp>::push_back() at 1 incorrect",
         std::vector< tk::real >{ 3.0, 4.0 }, p[1] );
  veceq( "<UnkEqComp>::push_back() at 2 incorrect",
         std::vector< tk::real >{ 5.0, 6.0 }, p[2] );
  veceq( "<UnkEqComp>::push_back() at 3 incorrect",
         std::vector< tk::real >{ 0.2, 0.3 }, p[3] );

  // tk::Data::push_back() unimplemented with EqCompUnk data layout
}

//! Test tk::Data::resize()
template<> template<>
void Data_object::test< 38 >() {
  set_test_name( "resize" );

  // Test with UnkEqComp data layout
  tk::Data< tk::UnkEqComp > p( 3, 2 );
  p(0,0,0) = 1.0;  p(0,1,0) = 2.0;
  p(1,0,0) = 3.0;  p(1,1,0) = 4.0;
  p(2,0,0) = 5.0;  p(2,1,0) = 6.0;

  // Enlarge by 4*2 (nprop=2 stays constant)
  p.resize( p.nunk()+4 );

  ensure_equals( "nunk after <UnkEqComp>::resize() incorrect", p.nunk(), 7 );
  ensure_equals( "nprop after <UnkEqComp>::resize() incorrect", p.nprop(), 2 );

  using unittest::veceq;

  veceq( "<UnkEqComp>::resize() at 0 incorrect",
         std::vector< tk::real >{ 1.0, 2.0 }, p[0] );
  veceq( "<UnkEqComp>::resize() at 1 incorrect",
         std::vector< tk::real >{ 3.0, 4.0 }, p[1] );
  veceq( "<UnkEqComp>::resize() at 2 incorrect",
         std::vector< tk::real >{ 5.0, 6.0 }, p[2] );

  veceq( "<UnkEqComp>::resize() at 3 incorrect",
         std::vector< tk::real >{ 0.0, 0.0 }, p[3] );
  veceq( "<UnkEqComp>::resize() at 4 incorrect",
         std::vector< tk::real >{ 0.0, 0.0 }, p[4] );
  veceq( "<UnkEqComp>::resize() at 5 incorrect",
         std::vector< tk::real >{ 0.0, 0.0 }, p[5] );
  veceq( "<UnkEqComp>::resize() at 6 incorrect",
         std::vector< tk::real >{ 0.0, 0.0 }, p[6] );

  // tk::Data::resize() unimplemented with EqCompUnk data layout
}

//! Test tk::Data::resize()
template<> template<>
void Data_object::test< 39 >() {
  set_test_name( "resize with non-default value" );

  // Test with UnkEqComp data layout
  tk::Data< tk::UnkEqComp > p( 3, 2 );
  p(0,0,0) = 1.0;  p(0,1,0) = 2.0;
  p(1,0,0) = 3.0;  p(1,1,0) = 4.0;
  p(2,0,0) = 5.0;  p(2,1,0) = 6.0;

  // Enlarge by 4*2 (nprop=2 stays constant)
  p.resize( p.nunk()+4, -2.13 );

  ensure_equals( "nunk after <UnkEqComp>::resize() incorrect", p.nunk(), 7 );
  ensure_equals( "nprop after <UnkEqComp>::resize() incorrect", p.nprop(), 2 );

  using unittest::veceq;

  veceq( "<UnkEqComp>::resize() at 0 incorrect",
         std::vector< tk::real >{ 1.0, 2.0 }, p[0] );
  veceq( "<UnkEqComp>::resize() at 1 incorrect",
         std::vector< tk::real >{ 3.0, 4.0 }, p[1] );
  veceq( "<UnkEqComp>::resize() at 2 incorrect",
         std::vector< tk::real >{ 5.0, 6.0 }, p[2] );

  veceq( "<UnkEqComp>::resize() at 3 incorrect",
         std::vector< tk::real >{ -2.13, -2.13 }, p[3] );
  veceq( "<UnkEqComp>::resize() at 4 incorrect",
         std::vector< tk::real >{ -2.13, -2.13 }, p[4] );
  veceq( "<UnkEqComp>::resize() at 5 incorrect",
         std::vector< tk::real >{ -2.13, -2.13 }, p[5] );
  veceq( "<UnkEqComp>::resize() at 6 incorrect",
         std::vector< tk::real >{ -2.13, -2.13 }, p[6] );

  // tk::Data::resize() unimplemented with EqCompUnk data layout
}

//! Test tk::Data::rm()
template<> template<>
void Data_object::test< 40 >() {
  set_test_name( "rm" );

  // Test with UnkEqComp data layout with a single component
  tk::Data< tk::UnkEqComp > p( 3, 1 );
  p(0,0,0) = 1.0;
  p(1,0,0) = 3.0;
  p(2,0,0) = 5.0;

  p.rm( { 0, 2 } );

  using unittest::veceq;

  ensure_equals( "nunk after <UnkEqComp>::rm(1prop) incorrect", p.nunk(), 1 );
  ensure_equals( "nprop after <UnkEqComp>::rm(1prop) incorrect", p.nprop(), 1 );
  veceq( "<UnkEqComp>::rm() incorrect", std::vector< tk::real >{ 3.0 }, p[0] );

  // Test with UnkEqComp data layout with a multiple components removing a
  // single unknown
  tk::Data< tk::UnkEqComp > q( 3, 2 );
  q(0,0,0) = 1.0;  q(0,1,0) = 2.0;
  q(1,0,0) = 3.0;  q(1,1,0) = 4.0;
  q(2,0,0) = 5.0;  q(2,1,0) = 6.0;

  q.rm( { 1 } );

  using unittest::veceq;

  ensure_equals( "nunk after <UnkEqComp>::rm(2prop) incorrect", q.nunk(), 2 );
  ensure_equals( "nprop after <UnkEqComp>::rm(2prop) incorrect", q.nprop(), 2 );
  veceq( "<UnkEqComp>::rm() at 0 incorrect",
         std::vector< tk::real >{ 1.0, 2.0 }, q[0] );
  veceq( "<UnkEqComp>::rm() at 1 incorrect",
         std::vector< tk::real >{ 5.0, 6.0 }, q[1] );

  // Test with UnkEqComp data layout with a multiple components removing
  // multiple unknowns
  tk::Data< tk::UnkEqComp > r( 3, 2 );
  r(0,0,0) = 1.0;  r(0,1,0) = 2.0;
  r(1,0,0) = 3.0;  r(1,1,0) = 4.0;
  r(2,0,0) = 5.0;  r(2,1,0) = 6.0;

  r.rm( { 0, 2 } );

  using unittest::veceq;

  ensure_equals( "nunk after <UnkEqComp>::rm(2prop) incorrect", r.nunk(), 1 );
  ensure_equals( "nprop after <UnkEqComp>::rm(2prop) incorrect", r.nprop(), 2 );
  veceq( "<UnkEqComp>::rm() at 0 incorrect",
         std::vector< tk::real >{ 3.0, 4.0 }, r[0] );
}

} // tut::

#endif // test_Data_h
