//******************************************************************************
/*!
  \file      src/RegTest/RegSuite.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 10:18:48 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Regression test suite class definition
  \details   Regression test suite class definition.
*/
//******************************************************************************

#include <string>
#include <iostream>

#include <boost/format.hpp>
#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "charm++.h"
#include "charm.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include "pup.h"
#include "pup_stl.h"

#include "Tags.h"
#include "RegSuite.h"
#include "Regression.h"
#include "regression.decl.h"
#include "regtest.decl.h"

extern CProxy_Main mainProxy;

using regtest::RegSuite;

RegSuite::RegSuite( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tag::verbose >() ? std::cout : std::clog ),
  m_ntest( 0 ),
  m_nrun( 0 ),
  m_ncomplete( 0 ),
  m_nfail( 0 ),
  m_nskip( 0 ),
  m_nwarn( 0 ),
  m_nexcp( 0 )
//******************************************************************************
// Constructor
//! \param[in] cmdline Data structure storing data from the command-line parser
//! \author J. Bakosi
//******************************************************************************
{
  // Get group name string passed in by -g
  const auto grp = cmdline.get< tag::group >();

  // If only select groups to be run, see if there is any that will run
  bool work = false;
  if (grp.empty())
    work = true;
  else
    work = true;        // for now there is always work

  // Quit if there is no work to be done
  if (!work) {

    m_print.note( "\nNo regression test groups to be executed because no test "
                  "group names match '" + grp + "'.\n" );
    mainProxy.finalize( false );

  } else {

    m_print.endpart();
    m_print.part( "Regression test suite" );
    m_print.reghead( "Regression tests computed", grp );

    // Fire up all tests using the Charm++ runtime system
    const std::vector< std::string > tests {
      "./charmrun +p1 Main/unittest -g Exc",
      "./charmrun +p2 Main/unittest -g MKL",
      "./charmrun +p3 Main/unittest -g MPI",
      "./charmrun +p4 Main/unittest"
    };
    for (const auto& t : tests)
      if (grp.empty()) {                        // consider all test groups
        spawntest( t );
      } else if (t.find(grp) != std::string::npos) {
        // spawn only the tests that match the string specified via -g string
        spawntest( t );
      }

  }
}

void
RegSuite::spawntest( const std::string& t )
//******************************************************************************
//  Fire up a regression test
//! \param[in] t Name of the test
//! \author J. Bakosi
//******************************************************************************
{
  // Increase number of tests to run
  ++m_ntest;

  // Asynchronously fire up a test
  CProxy_Regression< CProxy_RegSuite >::ckNew( thisProxy, t );
}

void
RegSuite::evaluate( std::vector< std::string > status )
//******************************************************************************
// Evaluate a unit test
//! \param[in] status Vector strings containing the test results.
//! \author J. Bakosi
//******************************************************************************
{
  // Increase number of tests run
  ++m_nrun;

  // Evaluate test
  evaluate( status, m_ncomplete, m_nwarn, m_nskip, m_nexcp, m_nfail );

  // Echo one-liner info on result of test
  m_print.test( m_ncomplete, m_nfail, status );

  // Wait for all tests to finish, then quit
  if (m_nrun == m_ntest)
  {
//     assess( m_print, "serial and Charm++", m_nfail, m_nwarn, m_nskip, m_nexcp,
//             m_ncomplete );
    // Quit
    mainProxy.finalize( true );
  }
}

void
RegSuite::evaluate( std::vector< std::string > status,
                    std::size_t& ncomplete,
                    std::size_t& nwarn,
                    std::size_t& nskip,
                    std::size_t& nexcp,
                    std::size_t& nfail )
//******************************************************************************
//  Evaluate a single regression test
//! \param[in] status Vector of strings containing the test results. See
//!   regtest::Regression constructor for the expected structure of status.
//! \param[inout] ncomplete Number of completed tests
//! \param[inout] nwarn Number of tests with a warning
//! \param[inout] nskip Number of skipped tests
//! \param[inout] nexcp Number of tests with an exception
//! \param[inout] nfail Number of failed tests
//! \author J. Bakosi
//******************************************************************************
{
  ++ncomplete;                      // count number of tests completed
//   if (status[2] == "3")             // count number of tests with a warning
//     ++nwarn;
//   else if (status[2] == "7")        // count number of skipped tests
//     ++nskip;
//   else if (status[2] == "2")        // count number of tests throwing
//     ++nexcp;
//   else if (status[2] != "0")        // count number of failed tests
//     ++nfail;
}

void
RegSuite::assess( const tk::Print& print,
                  std::size_t nfail,
                  std::size_t nwarn,
                  std::size_t nskip,
                  std::size_t nexcp,
                  std::size_t ncomplete )
//******************************************************************************
// Echo final assessment after the full regression test suite has finished
//! \param[in] print Pretty printer
//! \param[in] nfail Number of failed tests
//! \param[in] nwarn Number of tests with a warning
//! \param[in] nskip Number of skipped tests
//! \param[in] nexcp Number of tests with an exception
//! \param[in] ncomplete Number of completed tests
//! \author J. Bakosi
//******************************************************************************
{
//   if (!nfail && !nwarn && !nskip && !nexcp) {
//     print.note< tk::QUIET >
//       ( "All " + std::to_string(ncomplete) + " " + suite + " tests passed" );
//   } else {
//     std::string skip, warn, fail, excp;
//     if (nwarn) warn = "finished with a warning: " + std::to_string(nwarn);
//     if (nskip) skip = std::string(nwarn ? ", " : "") +
//                       "skipped: " + std::to_string(nskip);
//     if (nexcp) excp = std::string(nskip || nwarn ? ", " : "") +
//                       "threw exception: " + std::to_string(nexcp);
//     if (nfail) fail = std::string(nexcp || nskip || nwarn ?
//                       ", " : "") + "failed: " + std::to_string(nfail);
//     print.note< tk::QUIET >
//               ( "Of " + std::to_string(ncomplete) + " tests total: " +
//                 warn + skip + excp + fail );
//   }
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <regsuite.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
