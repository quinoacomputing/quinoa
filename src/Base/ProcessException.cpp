// *****************************************************************************
/*!
  \file      src/Base/ProcessException.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Process an exception
  \details   This file contains the implementation of processing an exception.
    Logically, it would make sense to put this into Exception.C, however,
    Exception.h is included by all who want to be able throw an exception (a
    potentially large number of files) and that would pull in the charm++.h as
    well as the mpi.h headers, which triggers a slew of compiler warnings. On
    the other hand, processing an exception is only done by the executables'
    main chares objects and a few select explicit main() routines that use MPI,
    which is a lot less than those which throw, so processing an exception is
    separated here.
*/
// *****************************************************************************

#include <cstdio>
#include <csignal>
#include <exception>
#include <cfenv>

#include "NoWarning/charm.hpp"
#include "NoWarning/mpi.hpp"

#include "Exception.hpp"
#include "ProcessException.hpp"

namespace tk {

void
signalHandler( int signum )
// *****************************************************************************
// Signal handler for multiple signals, SIGABRT, SIGSEGV, etc.
//! \param[in] signum Signal number
//! \see https://oroboro.com/stack-trace-on-crash
//! \details Signals caught:
//!    SIGABRT is generated when the program calls the abort() function, such as
//!            when an assert() triggers
//!    SIGSEGV is generated when the program makes an illegal memory access, such
//!            as reading unaligned memory, dereferencing a null pointer, reading
//!            memory out of bounds etc.
//!     SIGILL is generated when the program tries to execute a malformed
//!            instruction. This happens when the execution pointer starts reading
//!            non-program data, or when a pointer to a function is corrupted.
//!     SIGFPE is generated when executing an illegal floating point instruction,
//!            most commonly division by zero or floating point overflow.
// *****************************************************************************
{
  // associate each signal with a signal name string.
  const char* name = nullptr;
  switch( signum ) {
    case SIGABRT: name = "SIGABRT";  break;
    case SIGFPE:  name = "SIGFPE";  break;
    case SIGILL:  name = "SIGILL";   break;
    case SIGINT:  name = "SIGINT";  break;
    case SIGSEGV: name = "SIGSEGV";  break;
    case SIGTERM: name = "SIGTERM";   break;
  }

  // Echo what signal is caught
  if ( name )
    fprintf( stderr, "Caught signal %d (%s)\n", signum, name );
  else
    fprintf( stderr, "Caught signal %d\n", signum );

  // Get and display backtrace
  tk::Exception("Signal caught").handleException();

  // Tell the runtime system to exit with a nonzero exit code
  CkExit(1);
}

int
setSignalHandlers()
// *****************************************************************************
// Set signal handlers for multiple signals, SIGABRT, SIGSEGV, etc
//! \return Ignore, used for calling in a constructor's initializer list
// *****************************************************************************
{
  // override Charm++'s terminate handler
  std::set_terminate( [](){
    tk::Exception("Terminate was called").handleException();
    // Tell the runtime system to exit with a nonzero exit code
    CkExit(1);
  } );

  signal( SIGABRT, tk::signalHandler );
  signal( SIGFPE,  tk::signalHandler );
  signal( SIGILL,  tk::signalHandler );
  signal( SIGINT,  tk::signalHandler );
  signal( SIGSEGV, tk::signalHandler );
  signal( SIGTERM, tk::signalHandler );

  // This is commented at this time, because there is at least a single SIGFPE
  // that is generated from a place we have no control over, e.g., pthreads.
  // Will revisit in the future, because this should be enabled to detect and
  // terminate with a trace on floating point exceptions.
  //feenableexcept( FE_ALL_EXCEPT );

  return 0;
}

void
processExceptionCharm()
// *****************************************************************************
//  Process an exception from the Charm++ runtime system
//! \details See Josuttis, The C++ Standard Library - A Tutorial and Reference,
//!    2nd Edition, 2012.
// *****************************************************************************
{
  try {
    throw;      // rethrow exception to deal with it here
  }
  // Catch tk::Exception
  catch ( tk::Exception& qe ) {
    if (!CkMyPe()) qe.handleException();
  }
  // Catch std::exception and transform it into tk::Exception without
  // file:line:func information
  catch ( std::exception& se ) {
    tk::Exception qe( se.what() );
    if (!CkMyPe()) qe.handleException();
  }
  // Catch uncaught exception
  catch (...) {
    tk::Exception qe( "Non-standard exception" );
    if (!CkMyPe()) qe.handleException();
  }

  // Tell the runtime system to exit with a nonzero exit code
  CkExit(1);
}

} // tk::
