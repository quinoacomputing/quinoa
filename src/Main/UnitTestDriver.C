//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.C
  \author    J. Bakosi
  \date      Mon 11 Aug 2014 04:49:34 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     UnitTestDriver that drives the unit test suite
  \details   UnitTestDriver that drives the unit test suite
*/
//******************************************************************************

#include <UnitTestDriver.h>
#include <TUTSuite.h>
#include <unittest.decl.h>
#include <Handler.h>

using unittest::UnitTestDriver;

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

} // unittest::

UnitTestDriver::UnitTestDriver( const UnitTestPrint& print,
                                const ctr::CmdLine& cmdline ) :
  m_print( print )
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)
  try {

    m_print.endpart();
    m_print.part( "Factory" );

    // Instantiate and run unit test suite. We only support Template Unit Test
    // suites at this point, so no factory instantiation, simply fire up a
    // Charm++ chare TUTSuite, which fires up and evaluates all unit tests
    // included in Main/UnitTest.C. Store proxy handle in global-scope to make
    // it available to individual unit tests that fire up Charm++ chares so they
    // can call back to the test suite. Since this is called inside the main
    // chare constructor, the Charm++ runtime system distributes the handle
    // along with all other global-scope data.
    g_suiteProxy = CProxy_TUTSuite::ckNew( cmdline );

  } catch (...) { tk::processException(); }
}
