//******************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 01:02:06 PM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Template Unit Test suite class declaration
  \details   Template Unit Test suite class declaration. In principle there can
    be unit test suites other than this one which uses the Template Unit Test
    library.
*/
//******************************************************************************
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

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "tutsuite.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

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
      , { "Base/PUPUtil", 11 }
      , { "Base/Timer", 1 }
    };

    //! Fire up all tests in a test group
    void spawngrp( const std::string& g );
};

} // unittest::

#endif // TUTSuite_h
