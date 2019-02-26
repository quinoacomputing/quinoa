// *****************************************************************************
/*!
  \file      src/Main/UnitTestDriver.C
  \copyright 2012-2015, J. Bakosi, 2016-2019, Los Alamos National Security, LLC.
  \brief     Unit test driver
  \details   Unit test driver.
*/
// *****************************************************************************

#include "UnitTestPrint.h"
#include "UnitTestDriver.h"
#include "UnitTest/CmdLine/CmdLine.h"
#include "QuinoaConfig.h"
#include "NoWarning/tutsuite.decl.h"

using unittest::UnitTestDriver;

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

} // unittest::

UnitTestDriver::UnitTestDriver( const UnitTestPrint&,
                                const ctr::CmdLine& cmdline )
// *****************************************************************************
//  Constructor
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
// *****************************************************************************
{
  // Instantiate (on PE 0 ) and run unit test suite. We only support Template
  // Unit Test suites at this point, so no factory instantiation, simply fire up
  // a Charm++ chare TUTSuite, which fires up and evaluates all unit tests
  // included in Main/UnitTest.C. Store proxy handle in global-scope to make it
  // available to individual unit tests that fire up Charm++ chares so they can
  // call back to the test suite. Since this is called inside the main chare
  // constructor, the Charm++ runtime system distributes the handle along with
  // all other global-scope data.
  g_suiteProxy = CProxy_TUTSuite::ckNew( cmdline, 0 );
}
