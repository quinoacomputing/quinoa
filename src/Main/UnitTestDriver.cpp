// *****************************************************************************
/*!
  \file      src/Main/UnitTestDriver.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit test driver
  \details   Unit test driver.
*/
// *****************************************************************************

#include "UnitTestDriver.hpp"
#include "UnitTest/CmdLine/CmdLine.hpp"
#include "QuinoaConfig.hpp"
#include "NoWarning/tutsuite.decl.h"
#include "Writer.hpp"

using unittest::UnitTestDriver;

namespace unittest {

extern CProxy_TUTSuite g_suiteProxy;

} // unittest::

UnitTestDriver::UnitTestDriver( const ctr::CmdLine& cmdline, int )
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
