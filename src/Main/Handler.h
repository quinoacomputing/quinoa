//******************************************************************************
/*!
  \file      src/Main/Handler.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 03:40:45 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Handler functions
  \details   Handler functions.
             Note that including this header requires the upstream definition of
             mainProxy, i.e., it must be included after a Charm++-generated
             decl.h that declares the main charm module unique to an executable.
             This is not very elegant, but this way the code below is reused
             among different executables.
*/
//******************************************************************************
#ifndef Handler_h
#define Handler_h

#include <Exception.h>

extern CProxy_Main mainProxy;

namespace tk {

static void processException()
//******************************************************************************
// Process an exception
//! \details This function can be used to process an exception. The following
//! cases are handled: (1) tk::Exception, (2) std::exception converted to
//! tk::Exception with file,func,line information, (3) uncaught exceptions
//! converted to tk::Exception without file,func,line information.
//! Exception safety: no-throw guarantee: this function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {
    throw;      // rethrow exception to deal with it here
  }
    // Catch Quina::Exceptions
    catch ( Exception& qe ) {
      qe.handleException();
    }
    // Catch std::exception and transform it into Quinoa::Exception without
    // file:line:func information
    catch ( std::exception& se ) {
      Exception qe( se.what() );
      qe.handleException();
    }
    // Catch uncaught exception
    catch (...) {
      Exception qe( "Non-standard exception" );
      qe.handleException();
    }

  // Tell the runtime system to exit
  mainProxy.finalize();
}

} // tk

#endif // Handler_h
