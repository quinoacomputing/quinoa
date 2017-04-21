// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestMake_list.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/Make_list.h
  \details   Unit tests for Base/Make_list.h
*/
// *****************************************************************************
#ifndef test_Make_list_h
#define test_Make_list_h

#include "NoWarning/tut.h"

#include <boost/mpl/at.hpp>

#include "Make_list.h"

namespace tut {

//! All tests in group inherited from this base
struct Make_list_common {

  template< typename... Ts >
  struct Test {
    using list = typename tk::make_list< Ts... >::type;
  };

};

//! Test group shortcuts
using Make_list_group = test_group< Make_list_common, MAX_TESTS_IN_GROUP >;
using Make_list_object = Make_list_group::object;

//! Define test group
static Make_list_group Make_list( "Base/Make_list" );

//! Test definitions for group

//! \brief Test that tk::make_list can correctly convert a variadic template
//!   argument pack to boost::mpl::list
//! \author J. Bakosi
template<> template<>
void Make_list_object::test< 1 >() {
  set_test_name( "make_list" );

  using List = Test< int, tk::real, std::string, uint8_t >::list;
  boost::mpl::at< List, boost::mpl::int_<0> >::type first = 2;
  boost::mpl::at< List, boost::mpl::int_<1> >::type second = 1.23;
  boost::mpl::at< List, boost::mpl::int_<2> >::type third = "blah";
  boost::mpl::at< List, boost::mpl::int_<3> >::type fourth = 4U;

  ensure_equals( "variadic list item 0 correctly stored", first, 2 );
  ensure_equals( "variadic list item 1 correctly stored", second, 1.23 );
  ensure_equals( "variadic list item 2 correctly stored", third, "blah" );
  ensure_equals( "variadic list item 3 correctly stored", fourth, 4U );
}

} // tut::

#endif // test_Make_list_h
