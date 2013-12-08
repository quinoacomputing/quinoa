//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Sat 07 Dec 2013 04:59:12 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <RNGTestDriver.h>
#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/CmdLine/Parser.h>
#include <RNGTest/InputDeck/InputDeck.h>
#include <SmallCrush.h>
#include <Crush.h>
#include <BigCrush.h>
#include <MKLRNG.h>

using rngtest::RNGTestDriver;

RNGTestDriver::RNGTestDriver(int argc, char** argv, const tk::Print& print)
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

  // Create pretty printer for Quinoa
  m_print = std::unique_ptr< RNGTestPrint >( new RNGTestPrint(m_control) );
  m_paradigm = std::unique_ptr< tk::Paradigm >( new tk::Paradigm(print) );
  m_timer = std::unique_ptr< tk::Timer >( new tk::Timer );

  print.endpart();
  print.part("Factory");

  // Register random number generators
  quinoa::ctr::RNG rng;
  std::list< quinoa::ctr::RNGType > regRNG;
  initRNGFactory( m_RNGFactory, rng, regRNG, m_paradigm->nthreads(),
                  m_control->get< ctr::param, ctr::mklrng >() );
  print.list("Registered random number generators", rng, regRNG);

  // Bundle up essentials
  m_base = std::unique_ptr< Base >(
           new Base(*m_print, *m_paradigm, *m_control, *m_timer, m_RNGFactory) );

  //! Initialize factories
  initFactories( print );

  // Instantiate battery
  ctr::BatteryType b = m_control->get< ctr::selected, ctr::battery >();
  if (b != ctr::BatteryType::NO_BATTERY) {
    m_battery = std::unique_ptr< Battery >( m_batteryFactory[b]() );
  }

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
  add< SmallCrush >( m_batteryFactory, regBatt, battery,
                     ctr::BatteryType::SMALLCRUSH, *m_base );
  add< Crush >( m_batteryFactory, regBatt, battery,
                ctr::BatteryType::CRUSH, *m_base );
  add< BigCrush >( m_batteryFactory, regBatt, battery,
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

  if (!control.get<ctr::title>().empty()) {
    print.section("Title", control.get<ctr::title>());
  }

  if (m_battery) {
    print.Section<ctr::Battery, ctr::selected, ctr::battery>();
    m_battery->print();
    print.section("RNG(s) tested");
    print.Mklparams< quinoa::ctr::RNG,
                     quinoa::ctr::MKLUniformMethod,
                     quinoa::ctr::MKLGaussianMethod >
                   ( control.get<ctr::selected, ctr::rng>(),
                     control.get<ctr::param, ctr::mklrng>() );
    print.raw("\n");
  } else {
    print.note( "No RNG battery specified" );
  }
}

void
RNGTestDriver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
  //! Run battery (if any)
  if (m_battery) m_battery->run();
}
