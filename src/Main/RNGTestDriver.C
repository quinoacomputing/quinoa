//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Wed 09 Oct 2013 10:29:56 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <Macro.h>
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
// #ifdef MKL_CALLS
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
// #ifdef MKL_CALLS
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
// #ifdef MKL_CALLS
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

RNGTestDriver::RNGTestDriver(int argc, char** argv, Base& base) :
  m_base(base)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] base      Essentials
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate command line parser
  CmdLineParser cmdParser(argc, argv, base);

  // Instantiate input deck parser
  InputDeckParser idParser(base);

IGNORE(m_base);
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
