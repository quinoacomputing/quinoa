//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Mon 29 Jul 2013 10:32:14 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTestDriver that drives the random number generator test suite
  \details   RNGTestDriver that drives the random number generator test suite
*/
//******************************************************************************

#include <iostream>

#include <RNGTestDriver.h>
#include <Timer.h>
#include <Exception.h>

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
  std::cout << "RNGTestDriver::execute()" << std::endl;
}
