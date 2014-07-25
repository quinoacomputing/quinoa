//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 09:23:13 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <Factory.h>
#include <RNGTestDriver.h>
#include <RNGTest/CmdLine/Parser.h>
#include <RNGTest/InputDeck/Parser.h>
#include <TestU01Suite.h>
#include <rngtest.decl.h>
#include <Handler.h>

#ifdef HAS_MKL
#include <MKLRNG.h>
#endif

namespace rngtest {

extern ctr::InputDeck g_inputdeck;

} // rngtest::

using rngtest::RNGTestDriver;

RNGTestDriver::RNGTestDriver( const RNGTestPrint& print,
                              const ctr::CmdLine& cmdline ) :
  m_print( print )
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here
  try {

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

    m_print.part( "Factory" );

    // Register batteries
    BatteryFactory bf;
    using ctr::BatteryType;
    tk::recordCharmModel< Battery, TestU01Suite >( bf, BatteryType::SMALLCRUSH );
    tk::recordCharmModel< Battery, TestU01Suite >( bf, BatteryType::CRUSH );
    tk::recordCharmModel< Battery, TestU01Suite >( bf, BatteryType::BIGCRUSH );
    m_print.list< ctr::Battery >( "Registered batteries", bf );
    m_print.endpart();

    m_print.part( "Problem" );
    if ( !g_inputdeck.get< tag::title >().empty() )
      m_print.section( "Title", g_inputdeck.get< tag::title >() );

    // Instantiate and run battery
    const auto s = bf.find( g_inputdeck.get< tag::selected, tag::battery >() );
    if (s != end(bf)) {
      Battery( s->second() );
    } else Throw( "Battery not found in factory" );

  } catch (...) { tk::processException(); }
}
