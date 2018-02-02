// *****************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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

#include "UnitTestPrint.h"
#include "UnitTest/CmdLine/CmdLine.h"

#include "NoWarning/tutsuite.decl.h"

namespace unittest {

//! Template Unit Test unit test suite
class TUTSuite : public CBase_TUTSuite {

  public:
    //! Constructor
    explicit TUTSuite( const ctr::CmdLine& cmdline );

    //! Evaluate a unit test
    void evaluate( std::vector< std::string > status );

  private:
    UnitTestPrint m_print;      //!< Pretty printer
    std::size_t m_nrun;         //!< Number of tests ran (including dummies)
    std::size_t m_ngroup;       //!< Number of test groups
    std::size_t m_ncomplete;    //!< Number of completed tests
    std::size_t m_nfail;        //!< Number of failed tests
    std::size_t m_nskip;        //!< Number of skipped tests
    std::size_t m_nwarn;        //!< Number of tests with a warning
    std::size_t m_nexcp;        //!< Number of tests with an exception
    std::size_t m_nmigr;        //!< Number of Charm++ migration tests ran

    //! \brief Charm++ migration test group names and number of tests
    //! \details This map stores the names of test groups that define Charm++
    //!   migration tests and their associated number of Charm++ migration
    //!   tests. Every Charm++ migration test consists of two unit tests: one
    //!   for send and one for receive. Both triggers a TUT test, but the
    //!   receive side is created manually, i.e., without the awareness of the
    //!   TUT library. Unfortunately thus, there is no good way to count up
    //!   these additional tests, so they need to be explicitly maintained here.
    //!   To find out what tests spawn a new Charm++ chare, grep the src
    //!   directory for 'This test spawns a new Charm++ chare', which appears in
    //!   the comment before each Charm++ migration test name.
    const std::map< std::string, std::size_t > m_migrations {
        { "Base/Factory", 2 }
      , { "Base/PUPUtil", 14 }
      , { "Base/Timer", 1 }
      , { "Inciter/Scheme", 3 }
    };

    // Tests that must be run on PE 0
    // \details Some Charm++ tests must be run on PE 0 because they create
    // Charm++ chare arrays whose ckNew() must be called on PE 0.
    const std::unordered_set< std::string > m_fromPE0 {
        { "LoadBalance/LinearMap"}
      , { "LoadBalance/UnsMeshMap" }
      , { "Inciter/Scheme" }
    };

    //! Fire up all tests in a test group
    void spawngrp( const std::string& g );
};

} // unittest::

#endif // TUTSuite_h
