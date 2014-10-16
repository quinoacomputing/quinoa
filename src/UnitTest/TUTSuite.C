//******************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.C
  \author    J. Bakosi
  \date      Mon 04 Aug 2014 09:26:05 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Template Unit Test suite
  \details   Template Unit Test suite
*/
//******************************************************************************

#include <TUTSuite.h>
#include <TUTTest.h>
#include <unittest.decl.h>

extern CProxy_Main mainProxy;

namespace unittest {

extern tut::test_runner_singleton g_runner;

} // unittest::

using unittest::TUTSuite;

TUTSuite::TUTSuite( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tk::tag::verbose >() ? std::cout : std::clog ),
  m_maxTestsInGroup( 50 ),
  m_nrun( 0 ),
  m_ncomplete( 0 ),
  m_nfail( 0 ),
  m_nskip( 0 ),
  m_nwarn( 0 ),
  m_nexcp( 0 )
//******************************************************************************
// Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Output registered test groups
  const auto& groups = g_runner.get().list_groups();
  m_print.list( "Registered test groups", groups );
  m_ngroup = groups.size();

  m_print.endpart();
  m_print.part( "Problem" );
  m_print.unithead( "Unit tests computed", m_ngroup );

  // Fire up all tests in all groups
  for (const auto& g : groups)
    for (int t=1; t<=m_maxTestsInGroup; ++t)
      CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t );
}

void
TUTSuite::evaluate( std::vector< std::string > status )
//******************************************************************************
// Evaluate a unit test
//! \author  J. Bakosi
//******************************************************************************
{
  // Increase number tests run (including dummies)
  ++m_nrun;

  // Evaluate test
  if (status[2] != "8") {             // only care about non-dummy tests
    ++m_ncomplete;                    // count number of tests completed
    if (status[2] == "3")             // count number of tests with a warning
      ++m_nwarn;
    else if (status[2] == "7")        // count number of skipped tests
      ++m_nskip;
    else if (status[2] == "2")        // count number of tests throwing
      ++m_nexcp;
    else if (status[2] != "0")        // count number of failed tests
      ++m_nfail;
  }

  // Echo one-liner info on result of test
  m_print.test( m_ncomplete, m_nfail, status );

  if ( m_nrun == m_ngroup * m_maxTestsInGroup + 8 ) {
    // Echo final assessment
    if (!m_nfail && !m_nwarn && !m_nskip && !m_nexcp) {
      m_print.note< tk::QUIET >
                  ( "All " + std::to_string(m_ncomplete) + " tests passed" );
    } else {
      std::string skip, warn, fail, excp;
      if (m_nwarn) warn = "finished with a warning: " + std::to_string(m_nwarn);
      if (m_nskip) skip = std::string(m_nwarn ? ", " : "") +
                          "skipped: " + std::to_string(m_nskip);
      if (m_nexcp) excp = std::string(m_nskip || m_nwarn ? ", " : "") +
                          "threw exception: " + std::to_string(m_nexcp);
      if (m_nfail) fail = std::string(m_nexcp || m_nskip || m_nwarn ?
                          ", " : "") + "failed: " + std::to_string(m_nfail);
      m_print.note< tk::QUIET >
                  ( "Of " + std::to_string(m_ncomplete) + " tests total: " +
                    warn + skip + excp + fail );
    }
    // Quit
    mainProxy.finalize();
  }
}

#include <tutsuite.def.h>
