//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:15:29 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     UnitTestDriver that drives the unit test suite
  \details   UnitTestDriver that drives the unit test suite
*/
//******************************************************************************

#include <Print.h>
#include <Handler.h>
#include <Factory.h>
#include <UnitTestDriver.h>
#include <UnitTest/CmdLine/Parser.h>
#include <unittest.decl.h>

extern CProxy_Main mainProxy;

using unittest::UnitTestDriver;

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

  } catch (...) { tk::processException(); }
}

void
UnitTestDriver::execute()
//******************************************************************************
//  Run battery
//! \author J. Bakosi
//******************************************************************************
{
  // Quit
  mainProxy.finalize();  
}
