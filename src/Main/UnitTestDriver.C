//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.C
  \author    J. Bakosi
  \date      Thu 24 Jul 2014 08:13:23 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     UnitTestDriver that drives the unit test suite
  \details   UnitTestDriver that drives the unit test suite
*/
//******************************************************************************

#include <Print.h>
#include <Factory.h>
#include <UnitTestDriver.h>
#include <UnitTest/CmdLine/Parser.h>
#include <unittest.decl.h>
#include <Handler.h>

#include <tut/tut_reporter.hpp>

// Include unit test groups
#include <UnitTests/Base/flip_map.h>
#include <UnitTests/Base/make_list.h>
#include <UnitTests/Base/StrConvUtil.h>

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
  std::cout << "=============";
  tut::reporter reporter;
  tut::runner.get().set_callback( &reporter );
  tut::runner.get().run_tests();
  std::cout << "status: " << reporter.all_ok() << std::endl;
  std::cout << "=============\n";

  // Quit
  mainProxy.finalize();  
}
