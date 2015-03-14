//******************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.C
  \author    J. Bakosi
  \date      Thu 12 Mar 2015 10:17:42 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Template Unit Test suite class definition
  \details   Template Unit Test suite class definition. In principle there can
    be unit test suites other than this one which uses the Template Unit Test
    library.
*/
//******************************************************************************

#include <TUTSuite.h>
#include <Assessment.h>
#include <TUTTest.h>
#include <unittest.decl.h>

extern CProxy_Main mainProxy;

namespace unittest {

extern tut::test_runner_singleton g_runner;
extern int g_maxTestsInGroup;

} // unittest::

using unittest::TUTSuite;

TUTSuite::TUTSuite( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tag::verbose >() ? std::cout : std::clog ),
  m_nmpi( 0 ),
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
  // Output registered test groups
  const auto& groups = g_runner.get().list_groups();
  m_print.list( "Registered test groups", groups );
  m_ngroup = groups.size();

  m_print.endpart();
  m_print.part( "Serial and Charm++ unit test suite" );
  m_print.unithead( "Unit tests computed" );

  // Asynchronously fire up all tests in all groups using the Charm++ runtime
  // system
  for (const auto& g : groups)
    if (g.find("MPI") == std::string::npos)     // don't start MPI test groups
      for (int t=1; t<=g_maxTestsInGroup; ++t)
        CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t );
    else
      ++m_nmpi;
}

void
TUTSuite::evaluate( std::vector< std::string > status )
//******************************************************************************
// Evaluate a unit test
//! \param[in] status Vector strings containing the test results. See
//!   unittest::TUTTest constructor for the expected structure of status.
//! \author J. Bakosi
//******************************************************************************
{
  // Increase number tests run (including dummies)
  ++m_nrun;

  // Evaluate test
  unittest::evaluate( status, m_ncomplete, m_nwarn, m_nskip, m_nexcp, m_nfail );

  // Echo one-liner info on result of test
  m_print.test( m_ncomplete, m_nfail, status );

  // The magic number in the conditional is the number of Charm++ migration
  // tests. Every Charm++ migration test consists of two unit tests: one for
  // send and one for receive. Both triggers a TUT test, but the receive side is
  // created manually, i.e., without the awareness of the TUT library.
  // Unfortunately thus, there is no good way to count up these additional
  // tests.
  if ( m_nrun ==
       (m_ngroup-m_nmpi) * static_cast<std::size_t>(g_maxTestsInGroup) + 13 )
  {
    assess( m_print, "serial and Charm++", m_nfail, m_nwarn, m_nskip, m_nexcp,
            m_ncomplete );
    // Quit
    mainProxy.finalize();
  }
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <tutsuite.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
