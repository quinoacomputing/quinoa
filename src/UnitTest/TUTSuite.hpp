// *****************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Template Unit Test suite class declaration
  \details   Template Unit Test suite class declaration. In principle there can
    be unit test suites other than this one which uses the Template Unit Test
    library.
*/
// *****************************************************************************
#ifndef TUTSuite_h
#define TUTSuite_h

#include <vector>
#include <map>
#include <string>
#include <iosfwd>
#include <cstddef>
#include <cstring>

#include "UnitTestPrint.hpp"
#include "UnitTest/CmdLine/CmdLine.hpp"

#include "NoWarning/tutsuite.decl.h"
#include "NoWarning/mpirunner.decl.h"

namespace unittest {

//! Template Unit Test unit test suite
class TUTSuite : public CBase_TUTSuite {

  public:
    //! Constructor
    explicit TUTSuite( const ctr::CmdLine& cmdline );

    //! Evaluate a unit test
    void evaluate( std::vector< std::string > status );

  private:
    ctr::CmdLine m_cmdline;        //!< Command line user input
    //! MPI unit test runner nodegroup proxy
    CProxy_MPIRunner< CProxy_TUTSuite > m_mpirunner;
    std::size_t m_nrun;      //!< Number of tests ran (including dummies)
    std::size_t m_ngroup;    //!< Number of test groups
    std::size_t m_ncomplete; //!< Number of completed tests
    std::size_t m_nfail;     //!< Number of failed tests
    std::size_t m_nskip;     //!< Number of skipped tests
    std::size_t m_nwarn;     //!< Number of tests with a warning
    std::size_t m_nexcp;     //!< Number of tests with an exception
    std::size_t m_nspaw;     //!< Number of additionallly spawned tests ran

    //! \brief Charm++ test group names that spawn additional tests and number
    //!   of tests they spawn
    //! \details This map stores the names of test groups that define Charm++
    //!   tests, such as migration tests, that spawn multiple tests and their
    //!   associated number tests they additionally spawn. Every such test
    //!   consists of multiple unit tests: one for the host test and one or more
    //!   for receive(s). All such tests trigger/spawn additional TUT test(s).
    //!   The receive side of the spawed tests are created manually, i.e.,
    //!   without the awareness of the TUT library. Unfortunately, there is no
    //!   good way to count up these additionally spawned tests, so they need to
    //!   be explicitly maintained here. To find out what tests spawn a new
    //!   Charm++ chare, grep the src directory for 'This test spawns a new
    //!   Charm++ chare', which appears in the comment before each such
    //!   host test name.
    const std::map< std::string, std::size_t > m_nspawned {
        { "Base/Factory", 2 }
      , { "Base/PUPUtil", 14 }
      , { "Base/Timer", 1 }
      , { "Inciter/Scheme", 4 }
      , { "LinearSolver/ConjugateGradients", 2+4+2+4 }
    };

    // Tests that must be run on PE 0
    // \details Some Charm++ tests must be run on PE 0 because they create
    // Charm++ chare arrays whose ckNew() must be called on PE 0.
    const std::unordered_set< std::string > m_fromPE0 {
        { "LoadBalance/LinearMap"}
      , { "LoadBalance/UnsMeshMap" }
      , { "LinearSolver/ConjugateGradients" }
      , { "Inciter/Scheme" }
    };

    //! Fire up all tests in a test group
    void spawngrp( const std::string& g );

    //! Create pretty printer specialized to UnitTest
    //! \return Pretty printer
    UnitTestPrint printer() const {
      return UnitTestPrint(
        m_cmdline.logname( m_cmdline.get< tag::io, tag::screen >(),
                           m_cmdline.get< tag::io, tag::nrestart >() ),
        m_cmdline.get< tag::verbose >() ? std::cout : std::clog,
        std::ios_base::app );
    }
};

} // unittest::

#endif // TUTSuite_h
