// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Control/TestControl.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Control/Control
  \details   Unit tests for Control/Control
*/
// *****************************************************************************
#ifndef test_Control_h
#define test_Control_h

#include "NoWarning/tut.h"

#include "Control.h"

namespace tut {

//! Tags for unit tests
struct tag1 {};
struct tag2 {};
struct tag3 {};
struct tag4 {};
struct tag5 {};
struct tag6 {};
struct tag7 {};
struct tag8 {};
struct tag9 {};

//! All tests in group inherited from this base
struct Control_common {
  using map_value_tuple = tk::tuple::tagged_tuple< tag1, tk::real >;
  using tuple1 = tk::tuple::tagged_tuple<
                   tag1, std::string,
                   tag2, int,
                   tag5, std::vector< int >,
                   tag6, std::vector< std::vector< int > >,
                   tag9, std::vector< std::vector< std::vector< int > > >,
                   tag7, std::map< int, std::string >,
                   tag8, std::map< int, map_value_tuple > >;
  using tuple2 = tk::tuple::tagged_tuple< tag1, std::string,
                   tag2, int,
                   tag3, tuple1,
                   tag5, std::vector< int >,
                   tag6, std::vector< std::vector< int > >,
                   tag9, std::vector< std::vector< std::vector< int > > >,
                   tag7, std::map< int, std::string >,
                   tag8, std::map< int, map_value_tuple > >;
  using control = tk::Control<
                   tag1, std::string,
                   tag2, int,
                   tag3, tuple1,
                   tag4, tuple2,
                   tag5, std::vector< int >,
                   tag6, std::vector< std::vector< int > >,
                   tag9, std::vector< std::vector< std::vector< int > > >,
                   tag7, std::map< int, std::string >,
                   tag8, std::map< int, map_value_tuple > >;
};

//! Test group shortcuts
using Control_group = test_group< Control_common, MAX_TESTS_IN_GROUP >;
using Control_object = Control_group::object;

//! Define test group
static Control_group Control( "Control/Control" );

//! Test definitions for group

//! Test set(), const-ref get() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 1 >() {
  set_test_name( "set(), const-ref get() at three depths" );

  control c;
  c.set< tag1 >( "blah1" );
  c.set< tag3, tag1 >( "blah2" );
  c.set< tag4, tag3, tag1 >( "blah3" );

  ensure_equals( "get() 1st level", c.get< tag1 >(), "blah1" );
  ensure_equals( "get() 2nd level", c.get< tag3, tag1 >(), "blah2" );
  ensure_equals( "get() 3rd level", c.get< tag4, tag3, tag1 >(), "blah3" );
}

//! Test set(), rvalue get() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 2 >() {
  set_test_name( "set(), rvalue get() at three depths" );

  control c;
  c.get< tag1 >() = "blah1";
  c.get< tag3, tag1 >() = "blah2";
  c.get< tag4, tag3, tag1 >() = "blah3";

  ensure_equals( "get() 1st level", c.get< tag1 >(), "blah1" );
  ensure_equals( "get() 2nd level", c.get< tag3, tag1 >(), "blah2" );
  ensure_equals( "get() 3rd level", c.get< tag4, tag3, tag1 >(), "blah3" );
}

//! Test store() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 3 >() {
  set_test_name( "store() at three depths" );

  control c;
  c.store< tag2 >( "1" );
  c.store< tag3, tag2 >( "2" );
  c.store< tag4, tag3, tag2 >( "3" );

  ensure_equals( "get() 1st level", c.get< tag2 >(), 1 );
  ensure_equals( "get() 2nd level", c.get< tag3, tag2 >(), 2 );
  ensure_equals( "get() 3rd level", c.get< tag4, tag3, tag2 >(), 3 );
}

//! Test push_back() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 4 >() {
  set_test_name( "push_back() at three depths" );

  control c;
  c.push_back< tag5 >( 1 );
  c.push_back< tag5 >( 1 );
  c.push_back< tag3, tag5 >( 2 );
  c.push_back< tag3, tag5 >( 2 );
  c.push_back< tag4, tag3, tag5 >( 3 );
  c.push_back< tag4, tag3, tag5 >( 3 );

  ensure_equals( "vector size at 1st level", c.get< tag5 >().size(), 2 );
  ensure_equals( "vector size at 2nd level", c.get< tag3, tag5 >().size(), 2 );
  ensure_equals( "vector size at 3rd level",
                 c.get< tag4, tag3, tag5 >().size(), 2 );

  ensure( "vector elements correct at 1st level",
          c.get< tag5 >() == std::vector< int >{ 1, 1 } );
  ensure( "vector elements correct at 2nd level",
          c.get< tag3, tag5 >() == std::vector< int >{ 2, 2 } );
  ensure( "vector elements correct at 3rd level",
          c.get< tag4, tag3, tag5 >() == std::vector< int >{ 3, 3 } );
}

//! Test push_back_back() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 5 >() {
  set_test_name( "push_back_back() at three depths" );

  control c;
  c.push_back< tag6 >();             // create an outer vector element
  c.push_back_back< tag6 >( 1 );     // push to inner vector
  c.push_back_back< tag6 >( 1 );
  c.push_back< tag3, tag6 >();
  c.push_back_back< tag3, tag6 >( 2 );
  c.push_back_back< tag3, tag6 >( 2 );
  c.push_back< tag4, tag3, tag6 >();
  c.push_back_back< tag4, tag3, tag6 >( 3 );
  c.push_back_back< tag4, tag3, tag6 >( 3 );

  ensure_equals( "vector size at 1st level",
                 c.get< tag6 >().back().size(), 2 );
  ensure_equals( "vector size at 2nd level",
                 c.get< tag3, tag6 >().back().size(), 2 );
  ensure_equals( "vector size at 3rd level",
                 c.get< tag4, tag3, tag6 >().back().size(), 2 );

  ensure( "vector elements correct at 1st level",
          c.get< tag6 >().back() == std::vector< int >{ 1, 1 } );
  ensure( "vector elements correct at 2nd level",
          c.get< tag3, tag6 >().back() == std::vector< int >{ 2, 2 } );
  ensure( "vector elements correct at 3rd level",
          c.get< tag4, tag3, tag6 >().back() == std::vector< int >{ 3, 3 } );
}

//! Test store_back() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 6 >() {
  set_test_name( "store_back() at three depths" );

  control c;
  c.store_back< tag5 >( "1" );
  c.store_back< tag5 >( "1" );
  c.store_back< tag3, tag5 >( "2" );
  c.store_back< tag3, tag5 >( "2" );
  c.store_back< tag4, tag3, tag5 >( "3" );
  c.store_back< tag4, tag3, tag5 >( "3" );

  ensure_equals( "vector size at 1st level", c.get< tag5 >().size(), 2 );
  ensure_equals( "vector size at 2nd level", c.get< tag3, tag5 >().size(), 2 );
  ensure_equals( "vector size at 3rd level",
                 c.get< tag4, tag3, tag5 >().size(), 2 );

  ensure( "vector elements correct at 1st level",
          c.get< tag5 >() == std::vector< int >{ 1, 1 } );
  ensure( "vector elements correct at 2nd level",
          c.get< tag3, tag5 >() == std::vector< int >{ 2, 2 } );
  ensure( "vector elements correct at 3rd level",
          c.get< tag4, tag3, tag5 >() == std::vector< int >{ 3, 3 } );
}

//! Test store_back_back() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 7 >() {
  set_test_name( "store_back_back() at three depths" );

  control c;
  c.push_back< tag6 >();                // create an outer vector element
  c.store_back_back< tag6 >( "1" );     // push to inner vector
  c.store_back_back< tag6 >( "1" );
  c.push_back< tag3, tag6 >();
  c.store_back_back< tag3, tag6 >( "2" );
  c.store_back_back< tag3, tag6 >( "2" );
  c.push_back< tag4, tag3, tag6 >();
  c.store_back_back< tag4, tag3, tag6 >( "3" );
  c.store_back_back< tag4, tag3, tag6 >( "3" );

  ensure_equals( "vector size at 1st level",
                 c.get< tag6 >().back().size(), 2 );
  ensure_equals( "vector size at 2nd level",
                 c.get< tag3, tag6 >().back().size(), 2 );
  ensure_equals( "vector size at 3rd level",
                 c.get< tag4, tag3, tag6 >().back().size(), 2 );

  ensure( "vector elements correct at 1st level",
          c.get< tag6 >().back() == std::vector< int >{ 1, 1 } );
  ensure( "vector elements correct at 2nd level",
          c.get< tag3, tag6 >().back() == std::vector< int >{ 2, 2 } );
  ensure( "vector elements correct at 3rd level",
          c.get< tag4, tag3, tag6 >().back() == std::vector< int >{ 3, 3 } );
}

//! Test store_back_back_back() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 8 >() {
  set_test_name( "store_back_back_back() at three depths" );

  control c;
  c.push_back< tag9 >();                 // create an outer-outer vector element
  c.push_back_back< tag9 >();            // create an outer vector element
  c.store_back_back_back< tag9 >( "1" ); // push to inner vector
  c.store_back_back_back< tag9 >( "1" );
  c.push_back< tag3, tag9 >();
  c.push_back_back< tag3, tag9 >();
  c.store_back_back_back< tag3, tag9 >( "2" );
  c.store_back_back_back< tag3, tag9 >( "2" );
  c.push_back< tag4, tag3, tag9 >();
  c.push_back_back< tag4, tag3, tag9 >();
  c.store_back_back_back< tag4, tag3, tag9 >( "3" );
  c.store_back_back_back< tag4, tag3, tag9 >( "3" );

  ensure_equals( "vector size at 1st level",
                 c.get< tag9 >().back().back().size(), 2 );
  ensure_equals( "vector size at 2nd level",
                 c.get< tag3, tag9 >().back().back().size(), 2 );
  ensure_equals( "vector size at 3rd level",
                 c.get< tag4, tag3, tag9 >().back().back().size(), 2 );

  ensure( "vector elements correct at 1st level",
          c.get< tag9 >().back().back() == std::vector< int >{ 1, 1 } );
  ensure( "vector elements correct at 2nd level",
          c.get< tag3, tag9 >().back().back() == std::vector< int >{ 2, 2 } );
  ensure( "vector elements correct at 3rd level",
          c.get< tag4, tag3, tag9 >().back().back() == std::vector< int >{ 3, 3 } );
}

//! Test insert() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 9 >() {
  set_test_name( "insert() at three depths" );

  control c;
  c.insert< tag7 >( 1, "one" );
  c.insert< tag7 >( 2, "two" );
  c.insert< tag7 >( 10 );                    // use default ctor for mapped_type
  c.insert< tag3, tag7 >( 3, "three" );
  c.insert< tag3, tag7 >( 4, "four" );
  c.insert< tag3, tag7 >( 10 );              // use default ctor for mapped_type
  c.insert< tag4, tag3, tag7 >( 5, "five" );
  c.insert< tag4, tag3, tag7 >( 6, "six" );
  c.insert< tag4, tag3, tag7 >( 10 );        // use default ctor for mapped_type

  ensure_equals( "map size at 1st level", c.get< tag7 >().size(), 3 );
  ensure_equals( "map size at 2nd level", c.get< tag3, tag7 >().size(), 3 );
  ensure_equals( "map size at 3rd level",
                 c.get< tag4, tag3, tag7 >().size(), 3 );

  ensure( "map elements correct at 1st level",
          c.get< tag7 >() ==
            std::map< int, std::string >{ {1,"one"}, {2,"two"}, {10,""} } );
  ensure( "map elements correct at 2nd level",
          c.get< tag3, tag7 >() ==
            std::map< int, std::string >{ {3,"three"}, {4,"four"}, {10,""} } );
  ensure( "map elements correct at 3rd level",
          c.get< tag4, tag3, tag7 >() ==
            std::map< int, std::string >{ {5,"five"}, {6,"six"}, {10,""} } );
}

//! Test insert_field() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 10 >() {
  set_test_name( "insert_field() at three depths" );

  control c;
  c.insert_field< int, tag1, tag8 >( 1, "1.2" );
  c.insert_field< int, tag1, tag8 >( 2, "3.14" );
  c.insert_field< int, tag1, tag3, tag8 >( 3, "2.1" );
  c.insert_field< int, tag1, tag3, tag8 >( 4, "3e-3" );
  c.insert_field< int, tag1, tag4, tag3, tag8 >( 5, "-2.3" );
  c.insert_field< int, tag1, tag4, tag3, tag8 >( 6, "-10" );

  ensure_equals( "map size at 1st level", c.get< tag8 >().size(), 2 );
  ensure_equals( "map size at 2nd level", c.get< tag3, tag8 >().size(), 2 );
  ensure_equals( "map size at 3rd level",
                 c.get< tag4, tag3, tag8 >().size(), 2 );

  ensure( "map elements correct at 1st level",
          c.get< tag8 >() ==
            std::map< int, map_value_tuple >{ {1,{1.2}}, {2,{3.14}} } );
  ensure( "map elements correct at 2nd level",
          c.get< tag3, tag8 >() ==
            std::map< int, map_value_tuple >{ {3,{2.1}}, {4,{3e-3}} } );
  ensure( "map elements correct at 3rd level",
          c.get< tag4, tag3, tag8 >() ==
            std::map< int, map_value_tuple >{ {5,{-2.3}}, {6,{-10}} } );
}

//! Test insert_opt() at three levels
//! \author J. Bakosi
template<> template<>
void Control_object::test< 11 >() {
  set_test_name( "insert_opt() at three depths" );

  control c;
  c.insert_opt< int, tag1, tk::real, tag8 >( 1, 1.2 );
  c.insert_opt< int, tag1, tk::real, tag8 >( 2, 3.14 );
  c.insert_opt< int, tag1, tk::real, tag3, tag8 >( 3, 2.1 );
  c.insert_opt< int, tag1, tk::real, tag3, tag8 >( 4, 3e-3 );
  c.insert_opt< int, tag1, tk::real, tag4, tag3, tag8 >( 5, -2.3 );
  c.insert_opt< int, tag1, tk::real, tag4, tag3, tag8 >( 6, -10 );

  ensure_equals( "map size at 1st level", c.get< tag8 >().size(), 2 );
  ensure_equals( "map size at 2nd level", c.get< tag3, tag8 >().size(), 2 );
  ensure_equals( "map size at 3rd level",
                 c.get< tag4, tag3, tag8 >().size(), 2 );

  ensure( "map elements correct at 1st level",
          c.get< tag8 >() ==
            std::map< int, map_value_tuple >{ {1,{1.2}}, {2,{3.14}} } );
  ensure( "map elements correct at 2nd level",
          c.get< tag3, tag8 >() ==
            std::map< int, map_value_tuple >{ {3,{2.1}}, {4,{3e-3}} } );
  ensure( "map elements correct at 3rd level",
          c.get< tag4, tag3, tag8 >() ==
            std::map< int, map_value_tuple >{ {5,{-2.3}}, {6,{-10}} } );
}

//! Test T convert( str ) feeding garbage
//! \author J. Bakosi
template<> template<>
void Control_object::test< 12 >() {
  set_test_name( "T convert(str) feeding garbage" );

  control c;

  try {
    c.convert< int >( "a" );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
}

//! Test T convert( str ) feeding convertible value
//! \author J. Bakosi
template<> template<>
void Control_object::test< 13 >() {
  set_test_name( "T convert(str)" );

  control c;

  try {
    int a = c.convert< int >( "345" );
    ensure_equals( "conversion", a, 345 );
  }
  catch ( tk::Exception& ) {
    fail( "should not throw exception" );
  }
}

//! Test str convert(T) feeding convertible value
//! \author J. Bakosi
template<> template<>
void Control_object::test< 14 >() {
  set_test_name( "str convert(T)" );

  control c;

  try {
    std::string a = c.convert( 345 );
    ensure_equals( "conversion", a, "345" );
  }
  catch ( tk::Exception& ) {
    fail( "should not throw exception" );
  }
}

} // tut::

#endif // test_Control_h
