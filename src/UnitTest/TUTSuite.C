//******************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.C
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 01:09:16 PM MDT
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
  m_maxTestsInGroup( 5 ),
  m_ncomplete( 0 )
//******************************************************************************
// Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Output registered test groups
  const auto& groups = g_runner.get().list_groups();
  m_print.list( "Registered test groups", groups );

  // Assign custom Template Unit Test reporter
  Reporter reporter;
  g_runner.get().set_callback( &reporter );

  // Fire up all tests in all groups
  for (const auto& g : groups)
    for (int t=0; t<m_maxTestsInGroup; ++t)
      CProxy_TUTTest< CProxy_TUTSuite >::ckNew( thisProxy, g, t );
}

void
TUTSuite::evaluate()
//******************************************************************************
// Evaluate a unit test
//! \author  J. Bakosi
//******************************************************************************
{
std::cout << "eval: " << ++m_ncomplete << std::endl;
//   m_print.test( ++m_ncomplete, 0 );

  if ( m_ncomplete > 14 ) {
    //std::cout << "status: " << reporter.all_ok() << std::endl;
    // Quit
    mainProxy.finalize();  
  }
}

#include <tutsuite.def.h>
