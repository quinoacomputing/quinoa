// *****************************************************************************
/*!
  \file      tests/unit/Control/Options/TestMKLGammaMethod.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Control/Options/MKLGammaMethod
  \details   Unit tests for Control/Options/MKLGammaMethod
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Options/MKLGammaMethod.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct MKLGammaMethod_common {
  MKLGammaMethod_common() : m() {}
  const tk::ctr::MKLGammaMethod m;
};

//! Test group shortcuts
using MKLGammaMethod_group =
  test_group< MKLGammaMethod_common, MAX_TESTS_IN_GROUP >;
using MKLGammaMethod_object = MKLGammaMethod_group::object;

//! Define test group
static MKLGammaMethod_group MKLGammaMethod( "Control/Options/MKLGammaMethod" );

//! Test definitions for group

//! Test that member function param() finds MKL parameter for method type
template<> template<>
void MKLGammaMethod_object::test< 1 >() {
  set_test_name( "param() finds MKL param" );
  ensure( "cannot find parameter",
          m.param( tk::ctr::MKLGammaMethodType::GNORM ) ==
            (VSL_RNG_METHOD_GAMMA_GNORM) );
}

//! Test that member function param() throws in DEBUG mode if can't find param
template<> template<>
void MKLGammaMethod_object::test< 2 >() {
  set_test_name( "param() throws" );

  try {

    m.param( static_cast< tk::ctr::MKLGammaMethodType >( 234 ) );
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
template<> template<>
void MKLGammaMethod_object::test< 3 >() {
  set_test_name( "copy constructor" );

  tk::ctr::MKLGammaMethod c( m );
  std::vector< tk::ctr::MKLGammaMethod > v;
  v.push_back( c );
  ensure( "copy constructor used to push_back a MKLGammaMethod object to a "
          "std::vector",
          v[0].param( tk::ctr::MKLGammaMethodType::GNORM_ACCURATE ) ==
            (VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE) );
}

//! Test move constructor
template<> template<>
void MKLGammaMethod_object::test< 4 >() {
  set_test_name( "move constructor" );

  tk::ctr::MKLGammaMethod c( m );
  std::vector< tk::ctr::MKLGammaMethod > v;
  v.emplace_back( std::move(c) );
  ensure( "move constructor used to emplace_back a MKLGammaMethod object to "
          "a std::vector",
           v[0].param( tk::ctr::MKLGammaMethodType::GNORM ) ==
             (VSL_RNG_METHOD_GAMMA_GNORM) );
}

//! Test copy assignment
template<> template<>
void MKLGammaMethod_object::test< 5 >() {
  set_test_name( "copy assignment" );

  tk::ctr::MKLGammaMethod c;
  c = m;
  ensure( "find param of copy-assigned MKLGammaMethod",
          c.param( tk::ctr::MKLGammaMethodType::GNORM ) ==
             (VSL_RNG_METHOD_GAMMA_GNORM) );
}

//! Test move assignment
template<> template<>
void MKLGammaMethod_object::test< 6 >() {
  set_test_name( "move assignment" );

  tk::ctr::MKLGammaMethod c;
  c = std::move( m );
  ensure( "find param of move-assigned MKLGammaMethod",
          c.param( tk::ctr::MKLGammaMethodType::GNORM_ACCURATE ) ==
             (VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE) );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
