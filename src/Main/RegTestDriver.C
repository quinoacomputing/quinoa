//******************************************************************************
/*!
  \file      src/Main/RegTestDriver.C
  \author    J. Bakosi
  \date      Sun 31 May 2015 06:33:50 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Regression test harness driver
  \details   Regression test harness driver.
*/
//******************************************************************************

#include <boost/format.hpp>

#include "RegTestPrint.h"
#include "RegTestDriver.h"

using regtest::RegTestDriver;

RegTestDriver::RegTestDriver( const RegTestPrint& print,
                              const ctr::CmdLine& cmdline ) :
  m_print( print ),
  m_cmdline( cmdline )
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
}

void
RegTestDriver::execute() const
//******************************************************************************
//  Execute regression test harness driver
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate and run regression test suite
  CProxy_RegSuite::ckNew( m_cmdline );
}
