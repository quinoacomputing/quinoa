// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestVector.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/Vector.h
  \details   Unit tests for Base/Vector.h
*/
// *****************************************************************************
#ifndef test_Vector_h
#define test_Vector_h

#include <unistd.h>

#include "NoWarning/tut.h"

#include "Vector.h"

namespace tut {

//! All tests in group inherited from this base
struct Vector_common {
  double precision = 1.0e-15;    // required floating-point precision
};

//! Test group shortcuts
using Vector_group = test_group< Vector_common, MAX_TESTS_IN_GROUP >;
using Vector_object = Vector_group::object;

//! Define test group
static Vector_group Vector( "Base/Vector" );

//! Test definitions for group

//! Test cross product
//! \author J. Bakosi
template<> template<>
void Vector_object::test< 1 >() {
  set_test_name( "cross product" );

  std::array< tk::real, 3 > v1{{ 3.0, -3.0, 1.0 }},
                            v2{{ 4.0, 9.0, 2.0 }},
                            correct_result{{ -15.0, -2.0, 39.0 }};

  const auto result = tk::cross( v1, v2 );
  ensure_equals( "cross product incorrect",
                 result[0], correct_result[0], precision );
  ensure_equals( "cross product incorrect",
                 result[1], correct_result[1], precision );
  ensure_equals( "cross product incorrect",
                 result[2], correct_result[2], precision );
}

//! Test cross product divided by scalar
//! \author J. Bakosi
template<> template<>
void Vector_object::test< 2 >() {
  set_test_name( "cross product divided by scalar" );

  std::array< tk::real, 3 > v1{{ 3.0, -3.0, 1.0 }},
                            v2{{ 4.0, 9.0, 2.0 }},
                            correct_result{{ -7.5, -1.0, 19.5 }};

  const auto result = tk::crossdiv( v1, v2, 2.0 );
  ensure_equals( "cross product divided by scalar incorrect",
                 result[0], correct_result[0], precision );
  ensure_equals( "cross product divided by scalar incorrect",
                 result[1], correct_result[1], precision );
  ensure_equals( "cross product divided by scalar incorrect",
                 result[2], correct_result[2], precision );
}

//! Test dot product
//! \author J. Bakosi
template<> template<>
void Vector_object::test< 3 >() {
  set_test_name( "dot product" );

  std::array< tk::real, 3 > v1{{ 1.0, 2.0, 3.0 }}, v2{{ 4.0, -5.0, 6.0 }};
  tk::real correct_result = 12.0;

  const auto result = tk::dot( v1, v2 );
  ensure_equals( "dot product incorrect", result, correct_result, precision );
}

//! Test triple product
//! \author J. Bakosi
template<> template<>
void Vector_object::test< 4 >() {
  set_test_name( "triple product" );

  std::array< tk::real, 3 > v1{{ -1.0, 3.0, 3.0 }},
                            v2{{ -2.0, 3.0, 1.0 }},
                            v3{{  0.0, 4.0, 0.0 }};
  tk::real correct_result = -20.0;

  const auto result = tk::triple( v1, v2, v3 );
  ensure_equals( "triple product incorrect", result, correct_result,
                 precision );
}

} // tut::

#endif // test_Vector_h
