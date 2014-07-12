//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Sun 06 Jul 2014 08:32:39 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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

#ifdef HAS_MKL
#include <MKLRNG.h>
#endif

namespace rngtest {

extern ctr::InputDeck g_inputdeck;

} // rngtest::

using rngtest::RNGTestDriver;

RNGTestDriver::RNGTestDriver( int argc, char** argv, const RNGTestPrint& print )
 : m_print( print )
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here
  try {

    // Parse command line into cmdline
    ctr::CmdLine cmdline;
    CmdLineParser cmdParser( argc, argv, m_print, cmdline );

    // Parse input deck into g_inputdeck, transfer cmdline (no longer needed)
    InputDeckParser inputdeckParser( m_print, cmdline, g_inputdeck );

    m_print.endpart();

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

    m_print.part("Factory");

    // Register batteries
    BatteryFactory bf;
    tk::recordCharmModel< Battery, TestU01Suite >
                        ( bf, ctr::BatteryType::SMALLCRUSH );
    tk::recordCharmModel< Battery, TestU01Suite >
                        ( bf, ctr::BatteryType::CRUSH );
    tk::recordCharmModel< Battery, TestU01Suite >
                        ( bf, ctr::BatteryType::BIGCRUSH );
    m_print.list< ctr::Battery >( "Registered batteries", bf );

    //! Echo information on random number generator test suite to be created
    m_print.endpart();
    m_print.part("Problem");

    if ( !g_inputdeck.get< tag::title >().empty() )
      m_print.section("Title", g_inputdeck.get< tag::title >());

    // Instantiate and run battery
    const auto s = bf.find( g_inputdeck.get< tag::selected, tag::battery >() );
    if (s != end(bf)) {
      Battery( s->second() );
    } else Throw( "Battery not found in factory" );

  } catch (...) { tk::processException(); }
}
