//******************************************************************************
/*!
  \file      src/Base/Handler.h
  \author    J. Bakosi
  \date      Tue 17 Jun 2014 07:19:06 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Handler functions
  \details   Handler functions
*/
//******************************************************************************
#ifndef Handler_h
#define Handler_h

#include <Exception.h>

namespace tk {

//! Process an exception
void processException();

//! Quinoa's own new handler
void newHandler();

//! Quinoa's own terminate handler
void terminateHandler();

// Quinoa's own unexpected handler
void unexpectedHandler();

} // tk

#endif // Handler_h
