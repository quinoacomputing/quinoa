//******************************************************************************
/*!
  \file      src/Base/ProcessException.C
  \author    J. Bakosi
  \date      Tue 24 Mar 2015 12:23:40 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Process an exception
  \details   This file contains the implementation of processing an exception.
    Logically, it would make sense to put this into Exception.C, however,
    Exception.h is included by all who want to be able throw an exception (a
    potentially large number of files) and that would pull in the charm++.h
    header, which triggers a slew of compiler warnings. On the other hand,
    processing an exception is only done by executable's main chares objects (a
    lot less than those which throw), so processing an exception is separated
    here.
*/
//******************************************************************************
#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <charm++.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <Exception.h>
#include <ProcessException.h>

namespace tk {

//******************************************************************************
//  Process an exception
//! \details See Josuttis, The C++ Standard Library - A Tutorial and Reference,
//!    2nd Edition, 2012.
//! \author J. Bakosi
//! \author J. Bakosi
//******************************************************************************
void processException() {
  try {
    throw;      // rethrow exception to deal with it here
  }
    // Catch Quina::Exceptions
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
  // Tell the Charm++ runtime system to exit
  CkExit();
}

} // tk::
