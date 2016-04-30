//******************************************************************************
/*!
  \file      src/UnitTest/TUTTest.h
  \author    J. Bakosi
  \date      Sat 30 Apr 2016 02:40:31 PM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Template Unit Test unit test class declaration
  \details   Template Unit Test unit test class declaration.
*/
//******************************************************************************
#ifndef TUTTest_h
#define TUTTest_h

#include <string>
#include <iosfwd>

#include <tut/tut.hpp>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
  #pragma GCC diagnostic ignored "-Wreorder"
#endif

#include "tuttest.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace unittest {

extern tut::test_runner_singleton g_runner;

//! \brief Generic Charm++ chare class for unit tests utilizing the Template
//!    Unit Test library
template< class Proxy >
class TUTTest : public CBase_TUTTest< Proxy > {

  public:
    //! Constructor: run test then call back to host proxy to evaluate it
    //! \param[in] proxy Host proxy to call back to after test has been run
    //! \param[in] groupname Name of the group the test belongs to
    //! \param[in] t Test number on test group
    explicit TUTTest( Proxy& proxy, const std::string& groupname, int t ) {
      tut::test_result tr;
      g_runner.get().run_test( groupname, t, tr );
      proxy.evaluate( { tr.group, tr.name, std::to_string(tr.result),
                        tr.message, tr.exception_typeid } );
    }
};

} // unittest::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
  #pragma GCC diagnostic ignored "-Wreorder"
#endif

#define CK_TEMPLATES_ONLY
#include "tuttest.def.h"
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // TUTTest_h
