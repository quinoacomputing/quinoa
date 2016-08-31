// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestDataLayout.h
  \author    J. Bakosi
  \date      Wed 31 Aug 2016 12:43:18 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/DataLayout.h
  \details   Unit tests for Base/DataLayout.h
*/
// *****************************************************************************
#ifndef test_DataLayout_h
#define test_DataLayout_h

#include <limits>
#include <array>
#include <vector>

#include "NoWarning/tut.h"

#include <boost/mpl/at.hpp>

#include "DataLayout.h"

namespace tut {

//! All tests in group inherited from this base
struct DataLayout_common {

  const tk::real prec = std::numeric_limits< tk::real >::epsilon();

  // Ensure equality of all element of a vector of reals
  void veceq( const std::string& msg,
              const std::vector< tk::real >& a,
              std::vector< tk::real >&& b )
  {
    std::transform( a.begin(), a.end(), b.begin(), b.begin(),
                    [ &msg, this ]( tk::real s, tk::real& d ){
                      ensure_equals( msg, s, d, this->prec ); return true; } );
  }

  // Ensure equality of all element of a array of reals
  template< std::size_t N >
  void veceq( const std::string& msg,
              const std::array< tk::real, N >& a,
              std::array< tk::real, N >&& b )
  {
    std::transform( a.begin(), a.end(), b.begin(), b.begin(),
                    [ &msg, this ]( tk::real s, tk::real& d ){
                      ensure_equals( msg, s, d, this->prec ); return true; } );
  }
};

//! Test group shortcuts
using DataLayout_group = test_group< DataLayout_common, MAX_TESTS_IN_GROUP >;
using DataLayout_object = DataLayout_group::object;

//! Define test group
static DataLayout_group DataLayout( "Base/DataLayout" );

//! Test definitions for group

//! Test that tk::DataLayout's constructor creates
//!   correctly sized arrays
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 1 >() {
  set_test_name( "correct size" );

  // Test all template specializations and all callable ways
  ensure_equals( "<UnkEqComp>::nunk() returns 2",
                 tk::DataLayout< tk::UnkEqComp >( 2, 3 ).nunk(), 2 );
  ensure_equals( "<UnkEqComp>::nprop() returns 3",
                 tk::DataLayout< tk::UnkEqComp >( 2, 3 ).nprop(), 3 );

  ensure_equals( "<EqCompUnk>::nunk() returns 2",
                 tk::DataLayout< tk::EqCompUnk >( 2, 3 ).nunk(), 2 );
  ensure_equals( "<EqCompUnk>::nprop() returns 3",
                 tk::DataLayout< tk::EqCompUnk >( 2, 3 ).nprop(), 3 );
}

//! Test that tk::DataLayout's operator() returns the correct value
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 2 >() {
  set_test_name( "operator() returns correct value" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 5 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 5 );

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

//! Test that tk::DataLayout's operator() throws for out-of-bounds indices
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 3 >() {
  set_test_name( "operator() throws for out-of-bounds indices" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 5 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 5 );

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

//! Test that tk::DataLayout's var( cptr() ) returns the correct value
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 4 >() {
  set_test_name( "var(cptr()) returns correct value" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 5 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 5 );

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

//! Test that tk::DataLayout's var() throws for out-of-bounds indices
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 5 >() {
  set_test_name( "var(cptr()) throws for out-of-bounds indices" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 5 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 5 );

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

//! Test that tk::DataLayout's var( cptr() ) is equivalent to operator()
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 6 >() {
  set_test_name( "var(cptr()) == operator()" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 6 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 6 );

  pp( 1, 2, 3 ) = 0.23;
  pe( 1, 2, 3 ) = 0.32;

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::var(cptr) value incorrect",
                 pp.var( pp.cptr(2,3), 1 ), pp(1,2,3), prec );
  ensure_equals( "<EqCompUnk>::var(cptr) value incorrect",
                 pe.var( pe.cptr(2,3), 1 ), pe(1,2,3), prec );
}

//! Test that tk::DataLayout's layou() returns correct string
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 7 >() {
  set_test_name( "layou()" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 6 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 6 );

  // Test all template specializations
  ensure_equals( "<UnkEqComp>::layou() correct", pp.layout(), "unknown-major" );
  ensure_equals( "<EqCompUnk>::layou() correct", pe.layout(), "equation-major" );
}

//! Test that tk::DataLayout's extract() returns correct vector of unknowns
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 8 >() {
  set_test_name( "extract() vector of unknowns" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 3 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 3 );

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

//! Test that tk::DataLayout's extract() returns correct vector components
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 9 >() {
  set_test_name( "extract() vector of components" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 3 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 3 );

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

//! Test that tk::DataLayout's operator[] returns correct vector components
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 10 >() {
  set_test_name( "operator[] to access vector of components" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 3 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 3 );

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

//! Test that tk::DataLayout's extract() returns correct array of four reals
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 11 >() {
  set_test_name( "extract() array of four reals for A,B,C,D" );

  tk::DataLayout< tk::UnkEqComp > pp( 8, 2 );
  tk::DataLayout< tk::EqCompUnk > pe( 8, 2 );

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

//! Test that tk::DataLayout's extract() returns correct array of four reals
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 12 >() {
  set_test_name( "extract() array of four reals for N[4]" );

  tk::DataLayout< tk::UnkEqComp > pp( 8, 2 );
  tk::DataLayout< tk::EqCompUnk > pe( 8, 2 );

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

//! Test that tk::DataLayout's fill() correctly fills complete data array
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 13 >() {
  set_test_name( "fill() all with the same value" );

  tk::DataLayout< tk::UnkEqComp > pp( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > pe( 3, 2 );

  pp.fill( 0.0 );
  pe.fill( 0.1 );

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

//! Test that tk::DataLayout's fill() correctly fills vector of unknowns
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 14 >() {
  set_test_name( "fill() vector of unknowns with the same value" );

  tk::DataLayout< tk::UnkEqComp > pp( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > pe( 3, 2 );

  pp.fill( 0, 0, 1.5 );
  pp.fill( 1, 0, 2.5 );

  pe.fill( 0, 0, 0.5 );
  pe.fill( 1, 0, -0.5 );

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

//! \brief Test that tk::DataLayout's memory layout, i.e., stores data with
//!   correct strides
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 15 >() {
  set_test_name( "strides" );

  tk::DataLayout< tk::UnkEqComp > pp( 2, 5 );
  tk::DataLayout< tk::EqCompUnk > pe( 2, 5 );

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

//! Test tk::DataLayout's copy constructor
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 16 >() {
  set_test_name( "copy constructor" );

  tk::DataLayout< tk::UnkEqComp > p( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e( 3, 2 );

  p.fill( 0.1 );
  e.fill( 0.2 );

  std::vector< tk::DataLayout< tk::UnkEqComp > > v;
  v.push_back( p );
  std::vector< tk::DataLayout< tk::EqCompUnk > > w;
  w.push_back( e );

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

//! Test tk::DataLayout's copy assignment
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 17 >() {
  set_test_name( "copy assignment" );

  tk::DataLayout< tk::UnkEqComp > p( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e( 3, 2 );

  p.fill( 0.1 );
  e.fill( 0.2 );

  tk::DataLayout< tk::UnkEqComp > p1( 3, 2 );
  p1 = p;

  tk::DataLayout< tk::EqCompUnk > e1( 3, 2 );
  e1 = e;

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

//! Test tk::DataLayout's move constructor
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 18 >() {
  set_test_name( "move constructor" );

  tk::DataLayout< tk::UnkEqComp > p( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e( 3, 2 );

  p.fill( 0.1 );
  e.fill( 0.2 );

  std::vector< tk::DataLayout< tk::UnkEqComp > > v;
  v.emplace_back( p );
  std::vector< tk::DataLayout< tk::EqCompUnk > > w;
  w.emplace_back( e );

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

//! Test tk::DataLayout's move assignment
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 19 >() {
  set_test_name( "move assignment" );

  tk::DataLayout< tk::UnkEqComp > p( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e( 3, 2 );

  p.fill( 0.1 );
  e.fill( 0.2 );

  auto p1 = std::move( p );
  auto e1 = std::move( e );

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

//! Test tk::DataLayout's operator-=
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 20 >() {
  set_test_name( "operator-=" );

  tk::DataLayout< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  p1 -= p2;
  e1 -= e2;

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
}

//! Test tk::DataLayout's operator-
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 21 >() {
  set_test_name( "operator-" );

  tk::DataLayout< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  auto p = p1 - p2;
  auto e = e1 - e2;

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
}

//! Test tk::DataLayout's operator+=
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 22 >() {
  set_test_name( "operator+=" );

  tk::DataLayout< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  p1 += p2;
  e1 += e2;

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
}

//! Test tk::DataLayout's operator+
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 23 >() {
  set_test_name( "operator+" );

  tk::DataLayout< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  auto p = p1 + p2;
  auto e = e1 + e2;

  // Test all template specializations
  veceq( "<UnkEqComp>::operator+() at 0,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p.extract( 0, 0 ) );
  veceq( "<UnkEqComp>::operator+() at 1,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p.extract( 1, 0 ) );
  veceq( "<UnkEqComp>::operator+() at 0,1 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, p.extract( 0, 1 ) );

  veceq( "<EqCompUnk>::operator+() at 0,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e.extract( 0, 0 ) );
  veceq( "<EqCompUnk>::operator+() at 1,0 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e.extract( 1, 0 ) );
  veceq( "<EqCompUnk>::operator+() at 0,1 incorrect",
         std::vector< tk::real >{ 0.4, 0.4, 0.4 }, e.extract( 0, 1 ) );
}

//! Test tk::DataLayout's operator*=
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 24 >() {
  set_test_name( "operator*=" );

  tk::DataLayout< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  p1 *= p2;
  e1 *= e2;

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
}

//! Test tk::DataLayout's operator*
//! \author J. Bakosi
template<> template<>
void DataLayout_object::test< 25 >() {
  set_test_name( "operator*" );

  tk::DataLayout< tk::UnkEqComp > p1( 3, 2 ), p2( 3, 2 );
  tk::DataLayout< tk::EqCompUnk > e1( 3, 2 ), e2( 3, 2 );

  p1.fill( 0.1 );       p2.fill( 0.3 );
  e1.fill( 0.3 );       e2.fill( 0.1 );

  auto p = p1 * p2;
  auto e = e1 * e2;

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
}

} // tut::

#endif // test_DataLayout_h
