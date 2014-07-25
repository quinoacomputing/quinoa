//******************************************************************************
/*!
  \file      src/UnitTest/TUTTest.h
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 04:55:00 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Template Unit Test unit test
  \details   Template Unit Test unit test
*/
//******************************************************************************
#ifndef TUTTest_h
#define TUTTest_h

#include <tut/tut.hpp>
#include <tuttest.decl.h>

namespace unittest {

extern tut::test_runner_singleton g_runner;

//! Template Unit Test unit test
template< class Proxy >
class TUTTest : public CBase_TUTTest< Proxy > {

  public:
    //! Constructor: run test then evaluate it
    explicit TUTTest( Proxy& proxy, const std::string& groupname, int t ) {
      tut::test_result tr;
      g_runner.get().run_test( groupname, t, tr );
      proxy.evaluate( { tr.group, tr.name, std::to_string(tr.result),
                        tr.message, tr.exception_typeid } );
    }
};

} // unittest::

#define CK_TEMPLATES_ONLY
#include <tuttest.def.h>
#undef CK_TEMPLATES_ONLY

#endif // TUTTest_h
