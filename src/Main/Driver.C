//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Thu Aug 29 15:34:25 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>
#include <Timer.h>
#include <Exception.h>

using namespace quinoa;

Driver::Driver()
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
try :
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

Driver::~Driver() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  finalize();
}

void
Driver::finalize() noexcept
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
