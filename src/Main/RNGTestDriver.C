//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Sat 25 Jan 2014 03:14:14 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <make_unique.h>

#include <Factory.h>
#include <RNGTestDriver.h>
#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/CmdLine/Parser.h>
#include <RNGTest/InputDeck/InputDeck.h>
#include <SmallCrush.h>
#include <Crush.h>
#include <BigCrush.h>
#include <MKLRNG.h>

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
  CmdLineParser cmdParser(argc, argv, print, cmdline);

  // Parse input deck into m_control, transfer cmdline (no longer needed)
  InputDeckParser inputdeckParser(print, std::move(cmdline), m_control);

  // Create pretty printer
  m_print = tk::make_unique< RNGTestPrint >( m_control );
  m_paradigm = tk::make_unique< tk::Paradigm >( print );
  m_timer = tk::make_unique< tk::Timer >();

  print.endpart();
  print.part("Factory");

  // Register random number generators
  tk::ctr::RNG rng;
  std::list< tk::ctr::RNGType > regRNG;
  initRNGFactory( m_RNGFactory, rng, regRNG, m_paradigm->nthreads(),
                  m_control->get< tag::param, tk::tag::mklrng >(),
                  m_control->get< tag::param, tk::tag::rngsse >() );
  print.list("Registered random number generators", rng, regRNG);

  // Bundle up essentials
  m_base = tk::make_unique< Base >( *m_print, *m_paradigm, *m_control,
                                    *m_timer, m_RNGFactory );

  //! Initialize factories
  initFactories( print );

  // Instantiate battery
  ctr::BatteryType b = m_control->get< tag::selected, tag::battery >();
  if (b != ctr::BatteryType::NO_BATTERY) {
    m_battery = std::unique_ptr< Battery >( m_batteryFactory[b]() );
  }

  // Query number of tests to run
  if (m_battery) m_ntest = m_battery->ntest();

  //! Echo information on random number generator test suite to be created
  echo();
}

void
RNGTestDriver::initFactories(const tk::Print& print)
//******************************************************************************
//  Initialize factories
//! \author  J. Bakosi
//******************************************************************************
{
  // Register batteries
  ctr::Battery battery;
  std::list< ctr::BatteryType > regBatt;
  tk::regist< SmallCrush >( m_batteryFactory, regBatt, battery,
                            ctr::BatteryType::SMALLCRUSH, *m_base );
  tk::regist< Crush >( m_batteryFactory, regBatt, battery,
                       ctr::BatteryType::CRUSH, *m_base );
  tk::regist< BigCrush >( m_batteryFactory, regBatt, battery,
                          ctr::BatteryType::BIGCRUSH, *m_base );
  print.list("Registered batteries", battery, regBatt);
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
    m_battery->print();
    print.section("RNGs tested");
    print.MKLParams( control.get< tag::selected, tk::tag::rng >(),
                     control.get< tag::param, tk::tag::mklrng >() );
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
  if (m_battery && m_ntest) m_battery->run();
}
