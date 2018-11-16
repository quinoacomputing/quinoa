// *****************************************************************************
/*!
  \file      src/Base/ProcessException.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Process an exception definition
  \details   This file contains the definition of processing an exception.
    Logically, it would make sense to put this into Exception.C, however,
    Exception.h is included by all who want to be able throw an exception (a
    potentially large number of files) and that would pull in the charm++.h
    header, which triggers a slew of compiler warnings. On the other hand,
    processing an exception is only done by executable's main chares objects (a
    lot less than those which throw), so processing an exception is separated
    here.
*/
// *****************************************************************************
#ifndef ProcessException_h
#define ProcessException_h

namespace tk {

//! Signal handler for multiple signals, SIGABRT, SIGSEGV, etc.
void signalHandler( int signum );

//! Set signal handlers for multiple signals, SIGABRT, SIGSEGV, etc
int setSignalHandlers();

//! Process an exception from the Charm++ runtime system
void processExceptionCharm();

//! Process an exception from the MPI runtime system
void processExceptionMPI();

} // tk::

#endif // ProcessException_h
