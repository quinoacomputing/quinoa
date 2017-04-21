// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Control/Options/TestMKLBetaMethod.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Control/Options/MKLBetaMethod
  \details   Unit tests for Control/Options/MKLBetaMethod
*/
// *****************************************************************************
#ifndef test_MKLBetaMethod_h
#define test_MKLBetaMethod_h

#include "NoWarning/tut.h"

#include "Options/MKLBetaMethod.h"

namespace tut {

//! All tests in group inherited from this base
struct MKLBetaMethod_common {
  MKLBetaMethod_common() : m() {}
  const tk::ctr::MKLBetaMethod m;
};

//! Test group shortcuts
using MKLBetaMethod_group =
  test_group< MKLBetaMethod_common, MAX_TESTS_IN_GROUP >;
using MKLBetaMethod_object = MKLBetaMethod_group::object;

//! Define test group
static MKLBetaMethod_group MKLBetaMethod( "Control/Options/MKLBetaMethod" );

//! Test definitions for group

//! Test that member function param() finds MKL parameter for method type
//! \author J. Bakosi
template<> template<>
void MKLBetaMethod_object::test< 1 >() {
  set_test_name( "param() finds MKL param" );
  ensure( "cannot find parameter",
          m.param( tk::ctr::MKLBetaMethodType::CJA ) ==
            (VSL_RNG_METHOD_BETA_CJA) );
}

//! Test that member function param() throws in DEBUG mode if can't find param
//! \author J. Bakosi
template<> template<>
void MKLBetaMethod_object::test< 2 >() {
  set_test_name( "param() throws" );

  try {

    m.param( static_cast< tk::ctr::MKLBetaMethodType >( 234 ) );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif

  } catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Cannot find parameter" ) !=
              std::string::npos );
  }
}

//! Test copy constructor
//! \author J. Bakosi
template<> template<>
void MKLBetaMethod_object::test< 3 >() {
  set_test_name( "copy constructor" );

  tk::ctr::MKLBetaMethod c( m );
  std::vector< tk::ctr::MKLBetaMethod > v;
  v.push_back( c );
  ensure( "copy constructor used to push_back a MKLBetaMethod object to a "
          "std::vector",
          v[0].param( tk::ctr::MKLBetaMethodType::CJA_ACCURATE ) ==
            (VSL_RNG_METHOD_BETA_CJA_ACCURATE) );
}

//! Test move constructor
//! \author J. Bakosi
template<> template<>
void MKLBetaMethod_object::test< 4 >() {
  set_test_name( "move constructor" );

  tk::ctr::MKLBetaMethod c( m );
  std::vector< tk::ctr::MKLBetaMethod > v;
  v.emplace_back( std::move(c) );
  ensure( "move constructor used to emplace_back a MKLBetaMethod object to "
          "a std::vector",
           v[0].param( tk::ctr::MKLBetaMethodType::CJA ) ==
             (VSL_RNG_METHOD_BETA_CJA) );
}

//! Test copy assignment
//! \author J. Bakosi
template<> template<>
void MKLBetaMethod_object::test< 5 >() {
  set_test_name( "copy assignment" );

  tk::ctr::MKLBetaMethod c;
  c = m;
  ensure( "find param of copy-assigned MKLBetaMethod",
          c.param( tk::ctr::MKLBetaMethodType::CJA ) ==
             (VSL_RNG_METHOD_BETA_CJA) );
}

//! Test move assignment
//! \author J. Bakosi
template<> template<>
void MKLBetaMethod_object::test< 6 >() {
  set_test_name( "move assignment" );

  tk::ctr::MKLBetaMethod c;
  c = std::move( m );
  ensure( "find param of move-assigned MKLBetaMethod",
          c.param( tk::ctr::MKLBetaMethodType::CJA_ACCURATE ) ==
             (VSL_RNG_METHOD_BETA_CJA_ACCURATE) );
}

} // tut::

#endif // test_MKLBetaMethod_h
