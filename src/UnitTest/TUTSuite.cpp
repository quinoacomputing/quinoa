// *****************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Template Unit Test suite class definition
  \details   Template Unit Test suite class definition. In principle there can
    be unit test suites other than this one which uses the Template Unit Test
    library.
*/
// *****************************************************************************

#include <iostream>
#include <utility>
#include <map>

#include "NoWarning/format.hpp"

#include "NoWarning/tut_runner.hpp"

#include "Tags.hpp"
#include "TUTSuite.hpp"
#include "TUTTest.hpp"
#include "MPIRunner.hpp"
#include "Assessment.hpp"

#include "NoWarning/unittest.decl.h"

extern CProxy_Main mainProxy;

namespace unittest {

extern tut::test_runner_singleton g_runner;
extern int g_maxTestsInGroup;

} // unittest::

using unittest::TUTSuite;

TUTSuite::TUTSuite( const ctr::CmdLine& cmdline ) :
  m_cmdline( cmdline ),
  m_mpirunner(),
  m_nrun( 0 ),
  m_ngroup( 0 ),
  m_ncomplete( 0 ),
  m_nfail( 0 ),
  m_nskip( 0 ),
  m_nwarn( 0 ),
  m_nexcp( 0 ),
  m_nspaw( 0 )
// *****************************************************************************
// Constructor
//! \param[in] cmdline Data structure storing data from the command-line parser
// *****************************************************************************
{
  const auto& groups = g_runner.get().list_groups();

  { auto print = printer();
    print.part( "Factory" );
    // Output registered test groups
    print.list( "Registered test groups", groups );
  } // ensure print is destructed (cannot collide with that of evaluate)

  // Get group name string passed in by -g
  const auto grp = cmdline.get< tag::group >();

  // If only select groups to be run, see if there is any that will run
  bool work = false;
  if (grp.empty() ||
      std::any_of( groups.cbegin(), groups.cend(),
         [&grp]( const std::string& g )
         { return g.find(grp) != std::string::npos; } ))
    work = true;

  // Quit if there is no work to be done
  if (!work) {

    printer().note( "\nNo test groups to be executed because no test group "
                    "names match '" + grp + "'.\n" );
    mainProxy.finalize( true );

  } else {

    { auto print = printer();
      print.endpart();
      print.part( "Serial, Charm++, and MPI unit test suites" );
      print.unithead( "Unit tests computed", grp );
    } // ensure print is destructed (cannot collied with that of evaluate)

    // Create MPI unit test runner nodegroup
    m_mpirunner = CProxy_MPIRunner< CProxy_TUTSuite >::ckNew( thisProxy );

    // Fire up all tests in all groups using the Charm++ runtime system
    for (const auto& g : groups) {
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
// *****************************************************************************
{
  ++m_ngroup;         // increase number of test groups to run

  if (g.find("MPI") != std::string::npos)

    m_mpirunner.rungroup( g );

  else {

    // Add up number of additionally-spawned tests (this is so we know how many
    // to expect results from)
    const auto it = m_nspawned.find( g );
    if (it != m_nspawned.end()) m_nspaw += it->second;
  
    // Asynchronously fire up all tests in test group
    for (int t=1; t<=g_maxTestsInGroup; ++t) {
      auto i = m_fromPE0.find( g );
      if (i != end(m_fromPE0))
        CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t, 0 );
      else
        CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t );
    }
  }
}

void
TUTSuite::evaluate( std::vector< std::string > status )
// *****************************************************************************
// Evaluate a unit test
//! \param[in] status Vector strings containing the test results. See
//!   unittest::TUTTest constructor for the expected structure of status.
// *****************************************************************************
{
  // Increase number tests run (including dummies)
  ++m_nrun;

  // Evaluate test
  unittest::evaluate( status, m_ncomplete, m_nwarn, m_nskip, m_nexcp, m_nfail );

  auto print = printer();

  // Echo one-liner info on result of test
  print.test( m_ncomplete, m_nfail, status );

  // Wait for all tests to finish, then quit
  if (m_nrun == m_ngroup*static_cast<std::size_t>(g_maxTestsInGroup) + m_nspaw)
  {
    auto pass = assess(print, m_nfail, m_nwarn, m_nskip, m_nexcp, m_ncomplete);
    mainProxy.finalize( pass );
  }
}

#include "NoWarning/tutsuite.def.h"
