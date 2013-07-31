//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Wed Jul 31 09:50:37 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <iostream>
#include <sstream>

#include <mkl_vsl.h>

extern "C" {
  #include <unif01.h>
  #include <bbattery.h>
}

#include <RNGTestDriver.h>
#include <Timer.h>
#include <Exception.h>

VSLStreamStatePtr stream;
const int brng = VSL_BRNG_MCG59;
const unsigned int seed = 0;
const int method = VSL_RNG_METHOD_UNIFORM_STD;
const int n = 1;
const double a = 0.0;
const double b = 1.0;
double r;

void Quinoa::MKLErrChk(int vslerr)
//******************************************************************************
//  Special error handler for MKL
//! \param[in]  vslerr     Error code
//! \author  J. Bakosi
//******************************************************************************
{
  if (vslerr != VSL_STATUS_OK)
    try {

      std::stringstream s;
      s << "MKL VSL Error: code " << vslerr;
      Throw(ExceptType::FATAL, s.str());

    } catch (Exception&) {
        throw;
      }
      catch (std::exception& e) {
        Throw(ExceptType::FATAL, e.what());
      }
      catch (...) {
        Throw(ExceptType::UNCAUGHT, "non-standard exception");
      }
}

void initMKL()
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vslNewStream(&stream, brng, seed);
#else  // NDEBUG
  Quinoa::MKLErrChk(vslNewStream(&stream, brng, seed));
#endif // NDEBUG
#endif
}

void finalizeMKL()
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vslDeleteStream(&stream);
#else  // NDEBUG
  Quinoa::MKLErrChk(vslDeleteStream(&stream));
#endif // NDEBUG
#endif
}

double MKL_VSL()
//******************************************************************************
//  MKL VSL test
//! \author J. Bakosi
//******************************************************************************
{
#ifdef MKL_CALLS
#ifdef NDEBUG
  vdRngUniform(method, stream, n, &r, a, b);
#else  // NDEBUG
  Quinoa::MKLErrChk(vdRngUniform(method, stream, n, &r, a, b));
#endif // NDEBUG
#endif

  return r;
}

using namespace Quinoa;

RNGTestDriver::RNGTestDriver(int argc, char** argv)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \author J. Bakosi
//******************************************************************************
try :
  Driver(argc, argv),
  m_timer(nullptr)
{

  // Instantiate timer object
  m_timer = new(std::nothrow) Timer;
  ErrChk(m_timer != nullptr, ExceptType::FATAL,
         "Cannot allocate memory for timer object");

} // Roll back changes and rethrow on error
  catch (std::exception&) {
    finalize();
    throw;
  }
  // Catch uncaught exceptions
  catch (...) {
    finalize();
    Throw(ExceptType::UNCAUGHT, "Non-standard exception");
  }

RNGTestDriver::~RNGTestDriver() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  finalize();
}

void
RNGTestDriver::finalize() noexcept
//******************************************************************************
//  Finalize
//! \details Cleanup either at the end of business as usual or due to an
//!          exception. No-throw guarantee: this member function never throws
//!          exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  if (m_timer) { delete m_timer; m_timer = nullptr; }
}

void
RNGTestDriver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
  initMKL();

  const char* name = "MKL VSL test";

  unif01_Gen *gen;
  gen = unif01_CreateExternGen01(const_cast<char*>(name), MKL_VSL);
  bbattery_SmallCrush(gen);
  unif01_DeleteExternGen01(gen);

  finalizeMKL();
}
