//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 09:36:30 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <make_unique.h>

#include <Print.h>
#include <Handler.h>
#include <Factory.h>
#include <RNG.h>
#include <RNGTestDriver.h>
#include <RNGTest/CmdLine/Parser.h>
#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/InputDeck/InputDeck.h>
#include <TestU01Suite.h>
#include <SmallCrush.h>
#include <Crush.h>
#include <BigCrush.h>

#ifdef HAS_MKL
#include <MKLRNG.h>
#endif

namespace rngtest {

#ifdef HAS_MKL
extern tk::ctr::RNGMKLParameters g_rngmklparam;
#endif

extern tk::ctr::RNGSSEParameters g_rngsseparam;

extern std::vector< tk::ctr::RNGType > g_selectedrng;

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
  try {

    // Create simple pretty printer
    tk::Print print;

    // Parse command line into cmdline
    ctr::CmdLine cmdline;
    CmdLineParser cmdParser( argc, argv, print, cmdline );

    // Parse input deck into m_control, transfer cmdline (no longer needed)
    InputDeckParser inputdeckParser( print, cmdline, m_control );

    print.endpart();

    // Save RNG parameters to global-scope (needed for migrating the RNGs)
    #ifdef HAS_MKL
    g_rngmklparam = m_control.get< tag::param, tk::tag::rngmkl >();
    #endif
    g_rngsseparam = m_control.get< tag::param, tk::tag::rngsse >();

    // Save selected RNGs to global-scope (needed for migrating the RNGs)
    g_selectedrng = m_control.get< tag::selected, tk::tag::rng >();

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
    tk::recordModel< Battery, TestU01Suite< SmallCrush > >
                   ( bf, ctr::BatteryType::SMALLCRUSH );
    tk::recordModel< Battery, TestU01Suite< Crush > >
                   ( bf, ctr::BatteryType::CRUSH );
    tk::recordModel< Battery, TestU01Suite< BigCrush > >
                   ( bf, ctr::BatteryType::BIGCRUSH );
    print.list< ctr::Battery >( "Registered batteries", bf );

    //! Echo information on random number generator test suite to be created
    print.endpart();
    print.part("Problem");

    if ( !m_control.get< tag::title >().empty() )
      print.section("Title", m_control.get< tag::title >());

    // Instantiate and run battery
    const auto s = bf.find( m_control.get< tag::selected, tag::battery >() );
    if (s != end(bf)) {
      Battery battery( s->second() );

      // Query number of tests to run
      auto ntest = battery.ntest();

      if (ntest) {
        // Create RNGTest-specific pretty printer
        RNGTestPrint rngprint( m_control );

        rngprint.battery( ntest, battery.nstat() );
        battery.print( rngprint );
        rngprint.section("RNGs tested");
        #ifdef HAS_MKL
        rngprint.MKLParams( m_control.get< tag::selected, tk::tag::rng >(),
                            m_control.get< tag::param, tk::tag::rngmkl >() );
        #endif
        rngprint.RNGSSEParams( m_control.get< tag::selected, tk::tag::rng >(),
                               m_control.get< tag::param, tk::tag::rngsse >() );
        rngprint.raw('\n');
    
        //if (battery && m_ntest) //battery->run();
        //  CProxy_BatteryRunner::ckNew( *m_battery );

      } else Throw( tk::ExceptType::FATAL, "No tests in selected battery" );
    } else Throw( tk::ExceptType::FATAL, "Battery not found in factory" );

  } catch (...) { tk::processException(); }
}
