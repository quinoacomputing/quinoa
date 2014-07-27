//******************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.C
  \author    J. Bakosi
  \date      Sat 26 Jul 2014 06:26:43 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
  m_print( cmdline.get< tk::tag::verbose >() ? std::cout : tk::null ),
  m_maxTestsInGroup( 50 ),
  m_nrun( 0 ),
  m_ncomplete( 0 ),
  m_nfail( 0 )
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
  ++m_nrun;     // increase number tests run (including dummies)

  // Evaluate test
  if (status[2] != "8") {             // if not dummy
    ++m_ncomplete;                    // increase number of tests completed
    if (status[2] != "0") ++m_nfail;  // count number of failed tests
  }

  // Echo one-liner info on result of test
  m_print.test( m_ncomplete, m_nfail, status );

  if ( m_nrun == m_ngroup*m_maxTestsInGroup ) {
    // Echo final assessment
    if (m_nfail == 0) {
      m_print.note< tk::QUIET >
                  ( "All " + std::to_string(m_ncomplete) + " tests passed" );
    } else {
      m_print.note< tk::QUIET >
                  ( std::to_string(m_nfail) + " of " +
                    std::to_string(m_ncomplete) + " tests failed" );
    }
    // Quit
    mainProxy.finalize();
  }
}

#include <tutsuite.def.h>
