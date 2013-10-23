//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Tue 22 Oct 2013 08:52:40 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <RNGTestDriver.h>
#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/CmdLine/Parser.h>

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

using namespace rngtest;

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

  // Bundle up essentials
  m_base = std::unique_ptr< Base >(
             new Base(*m_print, *m_paradigm, *m_control, *m_timer) );

  print.endpart();
  print.part("Problem setup");
  print.section("Title", m_control->get<ctr::title>());
}

void
RNGTestDriver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
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
