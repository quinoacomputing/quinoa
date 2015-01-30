//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/make_list.h
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 03:11:35 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit tests for Base/make_list.h
  \details   Unit tests for Base/make_list.h
*/
//******************************************************************************
#ifndef test_make_list_h
#define test_make_list_h

#include <tut/tut.hpp>
#include <make_list.h>
#include <boost/mpl/at.hpp>

namespace tut {

//! All tests in group inherited from this base
struct make_list_common {

  template< typename... Ts >
  struct Test {
    using list = typename tk::make_list< Ts... >::type;
  };

};

//! Test group shortcuts
using make_list_group = test_group< make_list_common >;
using make_list_object = make_list_group::object;

//! Define test group
make_list_group make_list( "Base/make_list" );

//! Test definitions for group

//! Test that tk::make_list can correctly convert a variadic template argument
//! pack to boost::mpl::list
template<> template<>
void make_list_object::test< 1 >() {
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

#endif // test_make_list_h
