//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 09:31:49 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit test driver
  \details   Unit test driver.
*/
//******************************************************************************

#include <boost/format.hpp>

#include "UnitTestPrint.h"
#include "UnitTestDriver.h"
#include "UnitTest/CmdLine/CmdLine.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "tutsuite.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

using unittest::UnitTestDriver;

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

} // unittest::

UnitTestDriver::UnitTestDriver( const UnitTestPrint& print,
                                const ctr::CmdLine& cmdline ) :
  m_print( print )
//******************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)

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
}
