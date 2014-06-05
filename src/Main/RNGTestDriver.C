//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Thu 05 Jun 2014 07:32:50 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <make_unique.h>

#include <Macro.h>
#include <Factory.h>
#include <RNGTestDriver.h>
#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/CmdLine/Parser.h>
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

} //rngtest::

using rngtest::RNGTestDriver;

RNGTestDriver::RNGTestDriver(int argc, char** argv, const tk::Print& print) :
  m_ntest(0)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] print     Simple pretty printer
//! \author J. Bakosi
//******************************************************************************
{
  // Parse command line into cmdline
  std::unique_ptr< ctr::CmdLine > cmdline;
  CmdLineParser cmdParser( argc, argv, print, cmdline );

  // Parse input deck into m_control, transfer cmdline (no longer needed)
  InputDeckParser inputdeckParser( print, std::move(cmdline), m_control );

  // Create RNGTest-specific pretty printer
  m_print = tk::make_unique< RNGTestPrint >( m_control );

  // Create paradigm
  m_paradigm = tk::make_unique< tk::Paradigm >( print );

  // Create timer
  m_timer = tk::make_unique< tk::Timer >();

  print.endpart();
  print.part("Factory");

  // Save RNG parameters to global-scope (needed for migrating the RNGs)
  #ifdef HAS_MKL
  g_rngmklparam = m_control->get< tag::param, tk::tag::rngmkl >();
  #endif
  g_rngsseparam = m_control->get< tag::param, tk::tag::rngsse >();

  // Save selected RNGs to global-scope (needed for migrating the RNGs)
  g_selectedrng = m_control->get< tag::selected, tk::tag::rng >();

  // Bundle up essentials
  m_base = tk::make_unique< Base >( *m_print, *m_paradigm, *m_control,
                                    *m_timer, m_RNGFactory );

//   //! Initialize factories
//   initFactories( print );
// 
//   // Instantiate battery
//   ctr::BatteryType b = m_control->get< tag::selected, tag::battery >();
//   if (b != ctr::BatteryType::NO_BATTERY) {
//     m_battery = std::unique_ptr< Battery >( m_batteryFactory[b]() );
//   }
// 
//   // Query number of tests to run
//   if (m_battery) m_ntest = m_battery->ntest();
// 
//   //! Echo information on random number generator test suite to be created
//   echo();
}

void
RNGTestDriver::initFactories( const tk::Print& print )
//******************************************************************************
//  Initialize factories
//! \author  J. Bakosi
//******************************************************************************
{
  IGNORE(print);
  // Register batteries
//  ctr::Battery battery;
//  std::list< ctr::BatteryType > regBatt;
  tk::record< TestU01Suite< SmallCrush > >
            ( m_batteryFactory, ctr::BatteryType::SMALLCRUSH, *m_base );
  tk::record< TestU01Suite< Crush > >
            ( m_batteryFactory, ctr::BatteryType::CRUSH, *m_base );
  tk::record< TestU01Suite< BigCrush > >
            ( m_batteryFactory, ctr::BatteryType::BIGCRUSH, *m_base );
//  print.list("Registered batteries", battery, regBatt);
}

void
RNGTestDriver::echo()
//******************************************************************************
//  Echo information on random number generator test suite
//! \author J. Bakosi
//******************************************************************************
{
  const RNGTestPrint& print = m_base->print;
  const ctr::InputDeck& control = m_base->control;

  print.endpart();
  print.part("Problem");

  if ( !control.get< tag::title >().empty() ) {
    print.section("Title", control.get< tag::title >());
  }

  if (m_battery && m_ntest) {
    print.battery( m_ntest, m_battery->nstat() );
    m_battery->print( print );
    print.section("RNGs tested");
    #ifdef HAS_MKL
    print.MKLParams( control.get< tag::selected, tk::tag::rng >(),
                     control.get< tag::param, tk::tag::rngmkl >() );
    #endif
    print.RNGSSEParams( control.get< tag::selected, tk::tag::rng >(),
                        control.get< tag::param, tk::tag::rngsse >() );
    print.raw("\n");
  } else {
    print.note( "No RNG battery or no RNGs specified" );
  }
}

void
RNGTestDriver::execute() const
//******************************************************************************
//  Run battery
//! \author J. Bakosi
//******************************************************************************
{
//   if (m_battery && m_ntest) //m_battery->run();
//     CProxy_BatteryRunner::ckNew( *m_battery );
}
