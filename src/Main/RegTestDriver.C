//******************************************************************************
/*!
  \file      src/Main/RegTestDriver.C
  \author    J. Bakosi
  \date      Fri 20 Mar 2015 11:52:47 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Regression test harness driver
  \details   Regression test harness driver.
*/
//******************************************************************************

#include <RegTestDriver.h>

using regtest::RegTestDriver;

RegTestDriver::RegTestDriver( const RegTestPrint& print,
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
}
