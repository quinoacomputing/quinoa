//******************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.C
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 04:55:32 PM MDT
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
  m_ncomplete( 0 )
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
  ++m_ncomplete;
  m_print.test( status );

  if ( m_ncomplete == m_ngroup*m_maxTestsInGroup ) {
    //std::cout << "status: " << reporter.all_ok() << std::endl;
    // Quit
    mainProxy.finalize();  
  }
}

#include <tutsuite.def.h>
