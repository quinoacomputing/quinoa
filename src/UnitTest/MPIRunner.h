// *****************************************************************************
/*!
  \file      src/UnitTest/MPIRunner.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ nodegroup to run MPI unit tests
  \details   Charm++ nodegroup to run MPI unit tests.
*/
// *****************************************************************************
#ifndef MPIRunner_h
#define MPIRunner_h

#include <string>

#include "NoWarning/tut.h"
#include "NoWarning/mpirunner.decl.h"

namespace unittest {

extern tut::test_runner_singleton g_runner;
extern int g_maxTestsInGroup;

//! Generic Charm++ nodegroup chare class for running MPI unit tests
template< class Proxy >
class MPIRunner : public CBase_MPIRunner< Proxy > {

  public:
    //! Constrcutor: store host proxy
    MPIRunner( const Proxy& proxy ) : m_host(proxy) {}

    //! Fire up all tests in a test group
    void rungroup( const std::string& groupname ) {
      for (int t=1; t<=g_maxTestsInGroup; ++t) {      
        tut::test_result tr;
        // Calls MPI
        g_runner.get().run_test( groupname, t, tr );
        auto nd = CBase_MPIRunner< Proxy >::thisIndex;
        if (CkNodeFirst(nd) == 0) {  // only send one result back
          m_host.evaluate( { tr.group, tr.name, std::to_string(tr.result),
                             tr.message, tr.exception_typeid } );
        }
      }
    }

  private:
    Proxy m_host;       //!< Host proxy
};

} // unittest::

#define CK_TEMPLATES_ONLY
#include "NoWarning/mpirunner.def.h"
#undef CK_TEMPLATES_ONLY

#endif // MPIRunner_h
