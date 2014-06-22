//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Thu 19 Jun 2014 11:01:49 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <Print.h>
#include <Handler.h>
#include <Factory.h>
#include <RNGTestDriver.h>
#include <RNGTest/CmdLine/Parser.h>
#include <RNGTest/InputDeck/Parser.h>
#include <TestU01Suite.h>
#include <TestStack.h>

#ifdef HAS_MKL
#include <MKLRNG.h>
#endif

namespace rngtest {

extern ctr::InputDeck g_inputdeck;
extern TestStack g_testStack;

} // rngtest::

using rngtest::RNGTestDriver;

RNGTestDriver::RNGTestDriver( int argc, char** argv )
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here
  try {

    // Create simple pretty printer
    tk::Print print;

    // Parse command line into cmdline
    ctr::CmdLine cmdline;
    CmdLineParser cmdParser( argc, argv, print, cmdline );

    // Parse input deck into g_inputdeck, transfer cmdline (no longer needed)
    InputDeckParser inputdeckParser( print, cmdline, g_inputdeck );

    // Initialize statistical test stack
    g_testStack = TestStack();

    print.endpart();

  } catch (...) { tk::processException(); }
}

void
RNGTestDriver::execute()
//******************************************************************************
//  Run battery
//! \author J. Bakosi
//******************************************************************************
{
  try {

    // Create simple pretty printer
    tk::Print print;

    print.part("Factory");

    // Register batteries
    BatteryFactory bf;
    tk::recordCharmModel< Battery, TestU01Suite >
                        ( bf, ctr::BatteryType::SMALLCRUSH );
    tk::recordCharmModel< Battery, TestU01Suite >
                        ( bf, ctr::BatteryType::CRUSH );
    tk::recordCharmModel< Battery, TestU01Suite >
                        ( bf, ctr::BatteryType::BIGCRUSH );
    print.list< ctr::Battery >( "Registered batteries", bf );

    //! Echo information on random number generator test suite to be created
    print.endpart();
    print.part("Problem");

    if ( !g_inputdeck.get< tag::title >().empty() )
      print.section("Title", g_inputdeck.get< tag::title >());

    // Instantiate and run battery
    const auto s = bf.find( g_inputdeck.get< tag::selected, tag::battery >() );
    if (s != end(bf)) {
      Battery battery( s->second() );
      battery.run();
    } else Throw( tk::ExceptType::FATAL, "Battery not found in factory" );

  } catch (...) { tk::processException(); }
}
