//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Sat 04 May 2013 07:51:00 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class definition
  \details   Exception base class definition
*/
//******************************************************************************

#include <iostream>
#include <sstream>
#include <cstdio>

#include <execinfo.h>

#include <Driver.h>
#include <Exception.h>

using namespace Quinoa;

Exception::Exception(const ExceptType except,
                     const std::string& message,
                     int number) noexcept
//******************************************************************************
//  Constructor: generate error message
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
try :
  m_except(move(except)),
  m_message(move(message))
{

  stringstream s;
  s << m_message;
  if (number) s << number;
  m_message = s.str();

} // Catch std::exception
  catch (exception& se) {
    // Emit warning and continue
    cout << "RUNTIME ERROR in Exception constructor: " << se.what() << endl
         << "Continuing anyway..." << endl;
  }
  // Catch uncaught exceptions
  catch (...) {
    // Emit warning and continue
    cout << "UNKNOWN EXCEPTION in Exception constructor" << endl
         << "Continuing anyway..." << endl;
  }

void
Exception::trace() noexcept
//******************************************************************************
//  Echo call trace
//! \details No-throw guarantee: this member function never throws exceptions.
//!          For more information see the libc manual at
//!          http://www.gnu.org/software/libc/manual/html_node/Backtraces.html
//! \author J. Bakosi
//******************************************************************************
{
  // Obtain trace, requires stdio.h, execinfo.h
  void* callstack[128];
  int frames = backtrace(callstack, 128);
  char** strs = backtrace_symbols(callstack, frames);

  // Echo trace
  for (int i=0; i<frames; ++i) {
    printf("%s\n", strs[i]);
  }

  // allocated by backtrace_symbols()
  free(strs);
}

void
Exception::echo(const char* msg) noexcept
//******************************************************************************
//  Echo message
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  printf(">>> %s: %s\n>>> TRACE:\n", msg, what());
  trace();
}

ErrCode
Exception::handleException(Driver* const driver) noexcept
//******************************************************************************
//  Handle Exception: Print cumulative message and handle criticality
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  switch (m_except) {

    case WARNING:
      echo("WARNING");
      return NONFATAL;

    case CUMULATIVE:
      echo("CUMULATIVE ERROR");
      return NONFATAL;

    case ERROR:
      echo("ERROR");
      return NONFATAL;

    case FATAL:
      echo("FATAL ERROR");
      printf(">>> Attempting cleanup & graceful exit...\n");
      driver->finalize();
      return FATAL_ERROR;

    case RUNTIME:
      echo("RUNTIME ERROR");
      printf(">>> Attempting cleanup & graceful exit...\n");
      driver->finalize();
      return FATAL_ERROR;

    case UNCAUGHT:
      echo("UNKNOWN ERROR");
      printf(">>> Attempting cleanup & graceful exit...\n");
      driver->finalize();
      return FATAL_ERROR;

    default:
      return NONFATAL;
  }
}
