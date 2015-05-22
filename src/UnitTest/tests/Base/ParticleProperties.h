//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/ParticleProperties.h
  \author    J. Bakosi
  \date      Thu 21 May 2015 09:12:51 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit tests for Base/ParticleProperties.h
  \details   Unit tests for Base/ParticleProperties.h
*/
//******************************************************************************
#ifndef test_ParticleProperties_h
#define test_ParticleProperties_h

#include <tut/tut.hpp>
#include <ParticleProperties.h>
#include <boost/mpl/at.hpp>

namespace tut {

//! All tests in group inherited from this base
struct ParticleProperties_common {

  template< typename... Ts >
  struct Test {
    using list = typename tk::make_list< Ts... >::type;
  };

};

//! Test group shortcuts
using ParticleProperties_group =
  test_group< ParticleProperties_common, MAX_TESTS_IN_GROUP >;
using ParticleProperties_object = ParticleProperties_group::object;

//! Define test group
ParticleProperties_group ParticleProperties( "Base/ParticleProperties" );

//! Test definitions for group

//! \brief Test that tk::ParticleProperties' default constructor creates zero
//!   size array
//! \author J. Bakosi
template<> template<>
void ParticleProperties_object::test< 1 >() {
  set_test_name( "default ctor yields 0-size array" );

  // Test all template specializations
  ensure_equals( "<ParEqComp>::size() returns 0",
                 tk::ParticleProperties< tk::ParEqComp >().size(), 0 );
  ensure_equals( "<EqCompPar>::size() returns 0",
                 tk::ParticleProperties< tk::EqCompPar >().size(), 0 );
}

//! \brief Test that tk::ParticleProperties' never leaves its underlying array
//!   pointer uninitialized
//! \author J. Bakosi
template<> template<>
void ParticleProperties_object::test< 2 >() {
  set_test_name( "underlying pointer initialized" );

  // Test all template specializations
  ensure_not( "<ParEqComp>::ptr() not nullptr",
              tk::ParticleProperties< tk::ParEqComp >().ptr() == nullptr );
  ensure_not( "<EqCompPar>::ptr() not nullptr",
              tk::ParticleProperties< tk::EqCompPar >().ptr() == nullptr );
}

//! \brief Test that tk::ParticleProperties' non-default constructor creates
//!   correctly sized arrays
//! \author J. Bakosi
template<> template<>
void ParticleProperties_object::test< 3 >() {
  set_test_name( "ctor yields correct-size array" );

  // Test all template specializations and all callable ways
  ensure_equals( "<ParEqComp>::size() returns 0",
                 tk::ParticleProperties< tk::ParEqComp >( 2 ).size(), 0 );
  ensure_equals( "<ParEqComp>::npar() returns 2",
                 tk::ParticleProperties< tk::ParEqComp >( 2 ).npar(), 2 );

  ensure_equals( "<ParEqComp>::size() returns 0",
                 tk::ParticleProperties< tk::ParEqComp >( 2, 3 ).size(), 6 );
  ensure_equals( "<ParEqComp>::npar() returns 2",
                 tk::ParticleProperties< tk::ParEqComp >( 2, 3 ).npar(), 2 );

  ensure_equals( "<EqCompPar>::size() returns 0",
                 tk::ParticleProperties< tk::EqCompPar >( 2 ).size(), 0 );
  ensure_equals( "<EqCompPar>::npar() returns 2",
                 tk::ParticleProperties< tk::EqCompPar >( 2 ).npar(), 2 );

  ensure_equals( "<EqCompPar>::size() returns 0",
                 tk::ParticleProperties< tk::EqCompPar >( 2, 3 ).size(), 6 );
  ensure_equals( "<EqCompPar>::npar() returns 2",
                 tk::ParticleProperties< tk::EqCompPar >( 2, 3 ).npar(), 2 );
}

//! \brief Test that tk::ParticleProperties' operator() returns the correct
//!   value
//! \author J. Bakosi
template<> template<>
void ParticleProperties_object::test< 4 >() {
  set_test_name( "operator() returns correct val" );

  tk::ParticleProperties< tk::ParEqComp > pp( 2, 5 );
  tk::ParticleProperties< tk::EqCompPar > pe( 2, 5 );

  pp( 1, 2, 2 ) = 0.23;
  pe( 1, 2, 2 ) = 0.32;

  // Test all template specializations
  ensure_equals( "<ParEqComp>::() value correct", pp(1,2,2), 0.23 );
  ensure_equals( "<EqCompPar>::() value correct", pe(1,2,2), 0.32 );
}

//! \brief Test that tk::ParticleProperties' cvar( cptr() ) returns the correct
//!   value
//! \author J. Bakosi
template<> template<>
void ParticleProperties_object::test< 5 >() {
  set_test_name( "cvar(cptr()) returns correct val" );

  tk::ParticleProperties< tk::ParEqComp > pp( 2, 5 );
  tk::ParticleProperties< tk::EqCompPar > pe( 2, 5 );

  pp( 1, 2, 2 ) = 0.23;
  pe( 1, 2, 2 ) = 0.32;

  // Test all template specializations
  ensure_equals( "<ParEqComp>::cvar(cptr) value correct",
                 pp.cvar( pp.cptr(2,2), 1 ), 0.23 );
  ensure_equals( "<EqCompPar>::cvar(cptr) value correct",
                 pe.cvar( pe.cptr(2,2), 1 ), 0.32 );
}

//! \brief Test that tk::ParticleProperties' cvar( cptr() ) is equivalent to
//!    operator()
//! \author J. Bakosi
template<> template<>
void ParticleProperties_object::test< 6 >() {
  set_test_name( "cvar(cptr()) == operator()" );

  tk::ParticleProperties< tk::ParEqComp > pp( 2, 6 );
  tk::ParticleProperties< tk::EqCompPar > pe( 2, 6 );

  pp( 1, 2, 3 ) = 0.23;
  pe( 1, 2, 3 ) = 0.32;

  // Test all template specializations
  ensure_equals( "<ParEqComp>::cvar(cptr) value correct",
                 pp.cvar( pp.cptr(2,3), 1 ), pp(1,2,3) );
  ensure_equals( "<EqCompPar>::cvar(cptr) value correct",
                 pe.cvar( pe.cptr(2,3), 1 ), pe(1,2,3) );
}

//! \brief Test that tk::ParticleProperties' major() returns correct string
//! \author J. Bakosi
template<> template<>
void ParticleProperties_object::test< 7 >() {
  set_test_name( "major()" );

  tk::ParticleProperties< tk::ParEqComp > pp( 2, 6 );
  tk::ParticleProperties< tk::EqCompPar > pe( 2, 6 );

  // Test all template specializations
  ensure_equals( "<ParEqComp>::major() correct", pp.major(), "particle-major" );
  ensure_equals( "<EqCompPar>::major() correct", pe.major(), "equation-major" );
}
} // tut::

#endif // test_ParticleProperties_h
