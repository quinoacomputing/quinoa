// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Inciter/TestScheme.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Unit tests for Inciter/Scheme.h
  \details   Unit tests for Inciter/Scheme.h
*/
// *****************************************************************************
#ifndef test_Scheme_h
#define test_Scheme_h

#include <unistd.h>

#include "NoWarning/tut.h"

#include "Scheme.h"

namespace tut {

//! All tests in group inherited from this base
struct Scheme_common {};

//! Test group shortcuts
using Scheme_group = test_group< Scheme_common, MAX_TESTS_IN_GROUP >;
using Scheme_object = Scheme_group::object;

//! Define test group
static Scheme_group Scheme( "Inciter/Scheme" );

//! Test definitions for group

//! Test if ctor uses the correct underlying type
template<> template<>
void Scheme_object::test< 1 >() {
  set_test_name( "Scheme ctor which" );

  inciter::Scheme c( inciter::ctr::SchemeType::CG );
  ensure_equals( "Scheme underlying type", c.which(), 0 );
  inciter::Scheme d( inciter::ctr::SchemeType::DG );
  ensure_equals( "Scheme underlying type", d.which(), 1 );
}

//! Test if operator[] returns the correct underlying type
template<> template<>
void Scheme_object::test< 2 >() {
  set_test_name( "Scheme operator[] which" );

  inciter::Scheme c( inciter::ctr::SchemeType::CG );
  ensure_equals( "Scheme underlying element type", c.which_element(), 0 );
  inciter::Scheme d( inciter::ctr::SchemeType::DG );
  ensure_equals( "Scheme underlying element type", d.which_element(), 1 );
}

} // tut::

#endif // test_Scheme_h
