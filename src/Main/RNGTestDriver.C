//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Thu Nov 14 10:47:05 2013
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

// #include <iostream>
// #include <sstream>
// 
// #include <mkl_vsl.h>
// 
// extern "C" {
//   #include <unif01.h>
//   #include <bbattery.h>
// }
// 
// #include <Timer.h>
// #include <Exception.h>
// #include <MKLTest.h>
// 
// VSLStreamStatePtr stream;
// const int brng = VSL_BRNG_MCG59;
// const unsigned int seed = 0;
// const int method = VSL_RNG_METHOD_UNIFORM_STD;
// const int n = 1;
// const double a = 0.0;
// const double b = 1.0;
// double r;

// void
// MKLErrChk(int vslerr) 
// //******************************************************************************
// //  Special error handler for MKL
// //! \param[in]  vslerr     Error code
// //! \author  J. Bakosi
// //******************************************************************************
// {
//   if (vslerr != VSL_STATUS_OK)
//     try {
// 
//       std::stringstream s;
//       s << "MKL VSL Error: code " << vslerr;
//       Throw(Quinoa::ExceptType::FATAL, s.str());
// 
//     } catch (Exception&) {
//         throw;
//       }
//       catch (std::exception& e) {
//         Throw(Quinoa::ExceptType::FATAL, e.what());
//       }
//       catch (...) {
//         Throw(Quinoa::ExceptType::UNCAUGHT, "non-standard exception");
//       }
// }
// 
// void initMKL()
// {
// #ifdef HAS_MKL
// #ifdef NDEBUG
//   vslNewStream(&stream, brng, seed);
// #else  // NDEBUG
//   Quinoa::MKLErrChk(vslNewStream(&stream, brng, seed));
// #endif // NDEBUG
// #endif
// }
// 
// void finalizeMKL()
// {
// #ifdef HAS_MKL
// #ifdef NDEBUG
//   vslDeleteStream(&stream);
// #else  // NDEBUG
//   Quinoa::MKLErrChk(vslDeleteStream(&stream));
// #endif // NDEBUG
// #endif
// }
// 
// double MKL_VSL()
// //******************************************************************************
// //  MKL VSL test
// //! \author J. Bakosi
// //******************************************************************************
// {
// #ifdef HAS_MKL
// #ifdef NDEBUG
//   vdRngUniform(method, stream, n, &r, a, b);
// #else  // NDEBUG
//   Quinoa::MKLErrChk(vdRngUniform(method, stream, n, &r, a, b));
// #endif // NDEBUG
// #endif
// 
//   return r;
// }

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

  //! Echo information on random number generator test suite to be created
  echo();

  // Instantiate battery
  ctr::BatteryType b = m_control->get< ctr::selected, ctr::battery >();
  if (b != ctr::BatteryType::NO_BATTERY) {
    m_battery = std::unique_ptr< Battery >( m_batteryFactory[b]() );
  }

  // Echo 'unspecified' if not battery is unspecified
  if (!m_battery) {
    print.section("RNG battery", std::string("unspecified"));
  }
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
  print.section("Title", control.get<ctr::title>());
  print.Section<ctr::Battery, ctr::selected, ctr::battery>();

  print.Mklparams< quinoa::ctr::RNG,
                   quinoa::ctr::MKLUniformMethod,
                   quinoa::ctr::MKLGaussianMethod >
                 ( control.get<ctr::param, ctr::mklrng>() );
  print.endpart();
}

void
RNGTestDriver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
  //! Run battery (if any)
  if (m_battery) {
    m_battery->run();
  }

  //MKLTest mkl(control());

//   initMKL();
// 
//   const char* name = "MKL VSL test";
// 
//   unif01_Gen *gen;
//   gen = unif01_CreateExternGen01(const_cast<char*>(name), MKL_VSL);
//   bbattery_SmallCrush(gen);
//   unif01_DeleteExternGen01(gen);
// 
//   finalizeMKL();
}
