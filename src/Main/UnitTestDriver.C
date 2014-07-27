//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.C
  \author    J. Bakosi
  \date      Sat 26 Jul 2014 05:32:04 PM MDT
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

UnitTestDriver::UnitTestDriver( const UnitTestPrint& print,
                                const ctr::CmdLine& cmdline ) :
  m_print( print ),
  m_cmdline( cmdline )
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)
  try {

    m_print.endpart();

  } catch (...) { tk::processException(); }
}

void
UnitTestDriver::execute()
//******************************************************************************
//  Run unit test suite
//! \author J. Bakosi
//******************************************************************************
{
  m_print.part( "Factory" );

  // Instantiate and run unit test suite. We only support Template Unit Test
  // suites at this point, so no factory instantiation, simply fire up a Charm++
  // chare TUTSuite, which fires up and evaluates all unit tests included in
  // Main/UnitTest.C.
  CProxy_TUTSuite::ckNew( m_cmdline );
}
