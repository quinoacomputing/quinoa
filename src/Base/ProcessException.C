// *****************************************************************************
/*!
  \file      src/Base/ProcessException.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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

#include <exception>

#include "NoWarning/charm.h"
#include "NoWarning/mpi.h"

#include "Exception.h"
#include "ProcessException.h"

namespace tk {

void
processExceptionCharm()
// *****************************************************************************
//  Process an exception from the Charm++ runtime system
//! \details See Josuttis, The C++ Standard Library - A Tutorial and Reference,
//!    2nd Edition, 2012.
//! \author J. Bakosi
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

  // Tell the runtime system to exit with a message
  CkAbort( "Exception caught" );
}

void
processExceptionMPI()
// *****************************************************************************
//  Process an exception from the MPI runtime system
//! \details See Josuttis, The C++ Standard Library - A Tutorial and Reference,
//!    2nd Edition, 2012.
//! \author J. Bakosi
// *****************************************************************************
{
  int peid;

  #if defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wold-style-cast"
  #endif
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );  
  #if defined(__clang__)
    #pragma GCC diagnostic pop
  #endif

  try {
    throw;      // rethrow exception to deal with it here
  }
  // Catch tk::Exception
  catch ( tk::Exception& qe ) {
    if (peid == 0) qe.handleException();
  }
  // Catch std::exception and transform it into tk::Exception without
  // file:line:func information
  catch ( std::exception& se ) {
    tk::Exception qe( se.what() );
    if (peid == 0) qe.handleException();
  }
  // Catch uncaught exception
  catch (...) {
    tk::Exception qe( "Non-standard exception" );
    if (peid == 0) qe.handleException();
  }

  #if defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wold-style-cast"
  #endif
  // Tell the runtime system to exit with error code
  MPI_Abort( MPI_COMM_WORLD, tk::ErrCode::FAILURE );
  #if defined(__clang__)
    #pragma GCC diagnostic pop
  #endif
}

} // tk::
