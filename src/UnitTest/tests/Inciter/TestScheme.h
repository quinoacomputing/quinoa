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

//! Charm++ chare to test Pack/Unpack utilities during network migration
class Receiver : public CBase_Receiver {
  public:
    //! Constructor taking (and migrating) a Scheme
    explicit Receiver( const inciter::Scheme& s,
                       int expected,
                       std::string label )
    {
      // Create test result struct, assume test is ok
      tut::test_result tr( "Inciter/Scheme", 1,
                           "Charm:migrate Scheme(" + label + ") 2",
                           tut::test_result::result_type::ok );
      try {
        int actual = s.which();
        // Evaluate test
        ensure_equals( "Scheme underlying type different after migrated",
                       actual, expected );
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
};

//! Test definitions for group

//! Test if ctor uses the correct underlying type
template<> template<>
void Scheme_object::test< 1 >() {
  set_test_name( "ctor which" );

  inciter::Scheme c( inciter::ctr::SchemeType::CG );
  ensure_equals( "Underlying type", c.which(), 0 );
  inciter::Scheme d( inciter::ctr::SchemeType::DG );
  ensure_equals( "Underlying type", d.which(), 1 );
}

//! Test if operator[] returns the correct underlying type
template<> template<>
void Scheme_object::test< 2 >() {
  set_test_name( "operator[] which" );

  inciter::Scheme c( inciter::ctr::SchemeType::CG );
  ensure_equals( "Underlying element type", c.which_element(), 0 );
  inciter::Scheme d( inciter::ctr::SchemeType::DG );
  ensure_equals( "Underlying element type", d.which_element(), 1 );
}

//! Test Pack/Unpack of Scheme holding CProxy_CG
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void Scheme_object::test< 3 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate Scheme(CG) 1" );

  CProxy_Receiver::ckNew(
    inciter::Scheme( inciter::ctr::SchemeType::CG ), 0, "CG" );
}

//! Test Pack/Unpack of Scheme holding CProxy_DG
//! \details Every Charm++ migration test, such as this one, consists of two
//!   unit tests: one for send and one for receive. Both trigger a TUT test,
//!   but the receive side is created manually, i.e., without the awareness of
//!   the TUT library. Unfortunately thus, there is no good way to count up
//!   these additional tests, and thus if a test such as this is added to the
//!   suite this number must be updated in UnitTest/TUTSuite.h in
//!   unittest::TUTSuite::m_migrations.
template<> template<>
void Scheme_object::test< 4 >() {
  // This test spawns a new Charm++ chare. The "1" at the end of the test name
  // signals that this is only the first part of this test: the part up to
  // firing up an asynchronous Charm++ chare. The second part creates a new test
  // result, sending it back to the suite if successful. If that chare never
  // executes, the suite will hang waiting for that chare to call back.
  set_test_name( "Charm:migrate Scheme(DG) 1" );

  CProxy_Receiver::ckNew(
    inciter::Scheme( inciter::ctr::SchemeType::DG ), 1, "DG" );
}

} // tut::

#endif // test_Scheme_h
