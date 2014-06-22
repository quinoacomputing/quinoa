//******************************************************************************
/*!
  \file      src/Base/Handler.C
  \author    J. Bakosi
  \date      Tue 17 Jun 2014 07:18:52 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Handler functions
  \details   Handler functions
*/
//******************************************************************************

#include <Handler.h>
#include <rngtest.decl.h>

extern CProxy_Main mainProxy;

namespace tk {

void processException()
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
      Exception qe(ExceptType::RUNTIME, se.what());
      qe.handleException();
    }
    // Catch uncaught exception
    catch (...) {
      Exception qe(ExceptType::UNCAUGHT, "Non-standard exception");
      qe.handleException();
    }

  // Tell the runtime system to exit
  mainProxy.finalize();
}

void newHandler()
//******************************************************************************
// Quinoa's own new handler
//! \details This function is automatically called by allocation functions
//! whenever a memory allocation attempt fails. We simply throw an Exception
//! (derived from std::exception, that produces a trace. Exception safety:
//! no-throw guarantee: this function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {
    Throw(ExceptType::FATAL, "Cannot allocate memory");
  } catch (...) { processException(); }
}

void terminateHandler()
//******************************************************************************
// Quinoa's own terminate handler
//! \details This function is automatically called when the exception handling
//! process has to be abandoned for some reason. This happens when no catch
//! handler can be found for a thrown exception, or for some other exceptional
//! circumstance that makes impossible to continue the exception handling
//! process.We simply throw an Exception (derived from std::exception, that
//! produces a trace. Exception safety: no-throw guarantee: this function never
//! throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {
    Throw(ExceptType::FATAL, "Terminate called by uncaught exception");
  } catch (...) { processException(); }
}

void unexpectedHandler()
//******************************************************************************
// Quinoa's own unexpected handler
//! \details This function is automatically called when a function throws an
//! exception that is not in its dynamic-exception-specification (i.e., in its
//! throw specifier). The use of dynamic-exception-specifiers is deprecated
//! (since C++11). We simply throw an Exception (derived from std::exception,
//! that produces a trace. Exception safety: no-throw guarantee: this function
//! never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {
    Throw(ExceptType::FATAL,
          "Unexpected exception. "
          "Exception specifiers are deprecated since C++11.");
  } catch (...) { processException(); }
}

} // tk::
