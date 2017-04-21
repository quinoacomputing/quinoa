// *****************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Template Unit Test suite class definition
  \details   Template Unit Test suite class definition. In principle there can
    be unit test suites other than this one which uses the Template Unit Test
    library.
*/
// *****************************************************************************

#include <iostream>
#include <utility>
#include <map>

#include "NoWarning/format.h"

#include "NoWarning/tut_runner.h"

#include "Tags.h"
#include "TUTSuite.h"
#include "TUTTest.h"
#include "Assessment.h"

#include "NoWarning/unittest.decl.h"

extern CProxy_Main mainProxy;

namespace unittest {

extern tut::test_runner_singleton g_runner;
extern int g_maxTestsInGroup;

} // unittest::

using unittest::TUTSuite;

TUTSuite::TUTSuite( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tag::verbose >() ? std::cout : std::clog ),
  m_nrun( 0 ),
  m_ngroup( 0 ),
  m_ncomplete( 0 ),
  m_nfail( 0 ),
  m_nskip( 0 ),
  m_nwarn( 0 ),
  m_nexcp( 0 ),
  m_nmigr( 0 )
// *****************************************************************************
// Constructor
//! \param[in] cmdline Data structure storing data from the command-line parser
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.part( "Factory" );

  // Output registered test groups
  const auto& groups = g_runner.get().list_groups();
  m_print.list( "Registered test groups", groups );

  // Get group name string passed in by -g
  const auto grp = cmdline.get< tag::group >();

  // If only select groups to be run, see if there is any that will run
  bool work = false;
  if (grp.empty())
    work = true;
  else
    for (const auto& g : groups)
      if ( g.find("MPI") == std::string::npos &&   // don't consider MPI groups
           g.find(grp) != std::string::npos )
          work = true;

  // Quit if there is no work to be done
  if (!work) {

    m_print.note( "\nNo serial or Charm++ test groups to be executed because "
                  "no test group names match '" + grp + "'.\n" );
    mainProxy.finalize( false, true );

  } else {

    m_print.endpart();
    m_print.part( "Serial and Charm++ unit test suite" );
    m_print.unithead( "Unit tests computed", grp );

    // Fire up all tests in all groups using the Charm++ runtime system
    for (const auto& g : groups)
      if (g.find("MPI") == std::string::npos) {   // don't start MPI test groups
        if (grp.empty()) {                        // consider all test groups
          spawngrp( g );
        } else if (g.find(grp) != std::string::npos) {
          // spawn only the groups that match the string specified via -g string
          spawngrp( g );
        }
      }

  }
}

void
TUTSuite::spawngrp( const std::string& g )
// *****************************************************************************
//  Fire up all tests in a test group
//! \param[in] g Name of the test group
//! \author J. Bakosi
// *****************************************************************************
{
  ++m_ngroup;         // increase number of test groups to run

  // Add up number of Charm++ migration tests (this is so we know how many to
  // expect results from)
  const auto it = m_migrations.find( g );
  if (it != m_migrations.end()) m_nmigr += it->second;

  // Asynchronously fire up all tests in test group
  for (int t=1; t<=g_maxTestsInGroup; ++t) {
    auto i = m_fromPE0.find( g );
    if (i != end(m_fromPE0))
      CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t, 0 );
    else
      CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t );
  }
}

void
TUTSuite::evaluate( std::vector< std::string > status )
// *****************************************************************************
// Evaluate a unit test
//! \param[in] status Vector strings containing the test results. See
//!   unittest::TUTTest constructor for the expected structure of status.
//! \author J. Bakosi
// *****************************************************************************
{
  // Increase number tests run (including dummies)
  ++m_nrun;

  // Evaluate test
  unittest::evaluate( status, m_ncomplete, m_nwarn, m_nskip, m_nexcp, m_nfail );

  // Echo one-liner info on result of test
  m_print.test( m_ncomplete, m_nfail, status );

  // Wait for all tests to finish, then quit
  if (m_nrun == m_ngroup*static_cast<std::size_t>(g_maxTestsInGroup) + m_nmigr)
  {
    auto pass =
      assess( m_print, "serial and Charm++", m_nfail, m_nwarn, m_nskip, m_nexcp,
              m_ncomplete );
    // Quit
    mainProxy.finalize( true, pass );
  }
}

#include "NoWarning/tutsuite.def.h"
