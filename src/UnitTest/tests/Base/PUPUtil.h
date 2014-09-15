//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/PUPUtil.h
  \author    J. Bakosi
  \date      Sun 14 Sep 2014 10:31:17 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Unit tests for Base/PUPUtil.h
  \details   Unit tests for Base/PUPUtil.h
*/
//******************************************************************************
#ifndef test_PUPUtil_h
#define test_PUPUtil_h

#include <tut/tut.hpp>
#include <StrConvUtil.h>
#include <migrated.decl.h>
#include <tests/Base/MigratedTypes.h>

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

} // unittest::

namespace tut {

//! All tests in group inherited from this base
struct PUPUtil_common {};

//! Test group shortcuts
using PUPUtil_group = test_group< PUPUtil_common >;
using PUPUtil_object = PUPUtil_group::object;

//! Define test group
PUPUtil_group PUPUtil( "Base/PUPUtil" );

//! Test definitions for group

//! Charm++ chare to test Pack/Unpack utilities during network migration
struct Migrated : CBase_Migrated {

  //! Constructor taking (and migrating) a default strongly-typed enum
  explicit Migrated( charm::Enum_default e ) : m_enum_default(e) {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 1, 
                         "migrate default strongly-typed enum 2",
                         tut::test_result::result_type::ok );
    try {
      // Generate error message with expected and actual value in case if fail
      using underlying_type =
        typename std::underlying_type< charm::Enum_default >::type;
      std::string expected = std::to_string(
       static_cast< underlying_type >( charm::Enum_default::F1 ) );
      std::string actual = std::to_string(
       static_cast< underlying_type >( m_enum_default ) );
      // Evaluate test
      ensure( "strongly-typed enum different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_enum_default == charm::Enum_default::F1 );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluate(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a uint8_t strongly-typed enum
  explicit Migrated( charm::Enum_uint8_t e ) : m_enum_uint8_t(e) {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 2, 
                         "migrate uint_8t strongly-typed enum 2",
                         tut::test_result::result_type::ok );
    try {
      // Generate error message with expected and actual value in case if fail
      using underlying_type =
        typename std::underlying_type< charm::Enum_uint8_t >::type;
      std::string expected = std::to_string(
       static_cast< underlying_type >( charm::Enum_uint8_t::F1 ) );
      std::string actual = std::to_string(
       static_cast< underlying_type >( m_enum_uint8_t ) );
      // Evaluate test
      ensure( "strongly-typed enum different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_enum_uint8_t == charm::Enum_uint8_t::F1 );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluate(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a C-style enum
  explicit Migrated( charm::Enum_cstyle e ) : m_enum_cstyle(e) {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 3, 
                         "migrate C-style enum 2",
                         tut::test_result::result_type::ok );
    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = std::to_string( charm::Enum_cstyle::F1 );
      std::string actual = std::to_string( m_enum_cstyle );
      // Evaluate test
      ensure( "C-style enum different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_enum_cstyle == charm::Enum_cstyle::F1 );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluate(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::pair<int,double>
  explicit Migrated( charm::Pair p ) : m_pair(p) {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 4, 
                         "migrate std::pair<int,double> 2",
                         tut::test_result::result_type::ok );
    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "( 2, 3.14 )";
      std::string actual = "( " + std::to_string( m_pair.first ) +
                           ", " + std::to_string( m_pair.second ) + " )";
      // Evaluate test
      ensure( "std::pair different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_pair.first == 2 && m_pair.second == 3.14 );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluate(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::vector< std::string >
  explicit Migrated( charm::Vector v ) : m_vector(v) {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 5, 
                         "migrate std::vector< std::string > 2",
                         tut::test_result::result_type::ok );
    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "[ \"1\", \"blah\", \"boohoo\" ]";
      std::string actual = "[ " + m_vector[0] +
                           ", " + m_vector[1] +
                           ", " + m_vector[2] + " ]";
      // Evaluate test
      ensure( "std::pair different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_vector == charm::Vector{ "1", "blah", "boohoo" } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluate(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  //! Constructor taking (and migrating) a std::tuple
  explicit Migrated( charm::Tuple t ) : m_tuple(t) {
    // Create test result struct, assume test is ok
    tut::test_result tr( "Base/PUPUtil", 6, 
                         "migrate std::tuple 2",
                         tut::test_result::result_type::ok );
    try {
      // Generate error message with expected and actual value in case if fail
      std::string expected = "[ 2, 6.42, {\"woodoo\", \"boohoo\"}, "
                             "32, {{14, \"blah\"}, {15, \"blahblah\"}} ]";
      ensure_equals( "std::tuple different after migrated: vector size is "
                     "different", std::get< 2 >( m_tuple ).size(), 2UL );
      ensure_equals( "std::tuple different after migrated: map size is "
                     "different", std::get< 4 >( m_tuple ).size(), 2UL );
      using underlying_type_def =
        typename std::underlying_type< charm::Enum_default >::type;
      using underlying_type_uin =
        typename std::underlying_type< charm::Enum_uint8_t >::type;
      const auto& m = std::get< 4 >( m_tuple );
      const auto& i1 = m.find( charm::Enum_uint8_t::F1 );
      std::string m1, m2;
      if ( i1 != end(m) )
        m1 = "{" + std::to_string(
                     static_cast< underlying_type_uin >(
                       charm::Enum_uint8_t::F1 ) )
                 + ", \"" + i1->second + "\"}";
      else
        fail("std::tuple different after migrated: cannot find 1st key in map");
      const auto& i2 = m.find( charm::Enum_uint8_t::F2 );
      if ( i2 != end(m) )
        m2 = "{" + std::to_string(
                     static_cast< underlying_type_uin >(
                       charm::Enum_uint8_t::F2 ) )
                 + ", \"" + i2->second + "\"}";
      else
        fail("std::tuple different after migrated: cannot find 2nd key in map");
      std::string actual = "[ " + std::to_string( std::get< 0 >( m_tuple ) ) +
                           ", " + std::to_string( std::get< 1 >( m_tuple ) ) +
                        ", {\"" + std::get< 2 >( m_tuple )[ 0 ] +
                       "\", \"" + std::get< 2 >( m_tuple )[ 1 ] +
                        "\"}, " + std::to_string(
                                    static_cast< underlying_type_def >(
                                      std::get< 3 >( m_tuple ) ) ) +
                         ", {" + m1 + ", " + m2 + "} ]";
      // Evaluate test
      ensure( "std::tuple different after migrated: "
              "expected `" + expected + "` actual `" + actual + "`",
              m_tuple ==
                  charm::Tuple{ 2, 6.42, {"woodoo", "boohoo"},
                                charm::Enum_default::F1,
                                std::map< charm::Enum_uint8_t, std::string >{
                                  { charm::Enum_uint8_t::F1, "blah" },
                                  { charm::Enum_uint8_t::F2, "blahblah" } } } );
    } catch ( const failure& ex ) {
      tr.result = ex.result();
      tr.exception_typeid = ex.type();
      tr.message = ex.what();
    }
    // Send back a new test result, with tag "2", signaling the second part.
    unittest::g_suiteProxy.evaluate(
      { tr.group, tr.name, std::to_string(tr.result), tr.message,
        tr.exception_typeid } );
  }

  charm::Enum_default m_enum_default;
  charm::Enum_uint8_t m_enum_uint8_t;
  charm::Enum_cstyle m_enum_cstyle;
  charm::Pair m_pair;
  charm::Vector m_vector;
  charm::Tuple m_tuple;
};

//! Test Pack/Unpack of a default strongly-typed enum
template<> template<>
void PUPUtil_object::test< 1 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "migrate default strongly-typed enum 1" );

  CProxy_Migrated::ckNew( charm::Enum_default::F1 );
}

//! Test Pack/Unpack of a uint8_t strongly-typed enum
template<> template<>
void PUPUtil_object::test< 2 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "migrate uint8_t strongly-typed enum 1" );

  CProxy_Migrated::ckNew( charm::Enum_uint8_t::F1 );
}

//! Test Pack/Unpack of a C-style enum
template<> template<>
void PUPUtil_object::test< 3 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "migrate C-style enum 1" );

  CProxy_Migrated::ckNew( charm::Enum_cstyle::F1 );
}

//! Test Pack/Unpack of a std::pair<int,double>
template<> template<>
void PUPUtil_object::test< 4 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "migrate std::pair<int,double> 1" );

  CProxy_Migrated::ckNew( charm::Pair( 2, 3.14 ) );
}

//! Test Pack/Unpack of a std::vector< std::string >
template<> template<>
void PUPUtil_object::test< 5 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "migrate std::vector< std::string > 1" );

  CProxy_Migrated::ckNew( charm::Vector{"1", "blah", "boohoo"} );
}

//! Test Pack/Unpack of a std::tuple< int, double, std::vector< std::string >,
//! strongly-typed enum, std::map< strongly-typed(uint8_t) enum, std::string > >
template<> template<>
void PUPUtil_object::test< 6 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "migrate std::tuple 1" );

  CProxy_Migrated::ckNew(
    charm::Tuple{ 2, 6.42, {"woodoo", "boohoo"},
                  charm::Enum_default::F1,
                  std::map< charm::Enum_uint8_t, std::string >{
                    { charm::Enum_uint8_t::F1, "blah" },
                    { charm::Enum_uint8_t::F2, "blahblah" } } } );
}

} // tut::

#endif // test_PUPUtil_h
