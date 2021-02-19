// *****************************************************************************
/*!
  \file      src/UnitTest/MPIRunner.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ nodegroup to run MPI unit tests
  \details   Charm++ nodegroup to run MPI unit tests.
*/
// *****************************************************************************
#ifndef MPIRunner_h
#define MPIRunner_h

#include <string>

#include "NoWarning/tut.hpp"
#include "NoWarning/mpirunner.decl.h"

namespace unittest {

extern tut::test_runner_singleton g_runner;
extern int g_maxTestsInGroup;

//! Generic Charm++ nodegroup chare class for running MPI unit tests
template< class Proxy >
class MPIRunner : public CBase_MPIRunner< Proxy > {

  public:
    //! Constrcutor: store host proxy
    explicit MPIRunner( const Proxy& proxy ) : m_host(proxy) {}

    //! Fire up all tests in a test group
    void rungroup( const std::string& groupname ) {
      for (int t=1; t<=g_maxTestsInGroup; ++t) {      
        tut::test_result tr;
        auto nd = CBase_MPIRunner< Proxy >::thisIndex;
        // Invoke those MPI tests whose group name contains "MPISingle" from a
        // single MPI rank only
        if (groupname.find("MPISingle") != std::string::npos) {
          if (CkNodeFirst(nd) == 0) g_runner.get().run_test( groupname, t, tr );
        } else {
          g_runner.get().run_test( groupname, t, tr );
        }
        // Send result (only one)
        if (CkNodeFirst(nd) == 0) {
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
