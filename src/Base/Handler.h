//******************************************************************************
/*!
  \file      src/Base/Handler.h
  \author    J. Bakosi
  \date      Thu Sep 19 13:22:54 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Handler functions
  \details   Handler functions
*/
//******************************************************************************
#ifndef Handler_h
#define Handler_h

namespace quinoa {

//! Process an exception
void processException [[noreturn]] () noexcept;

//! Quinoa's own new handler
void newHandler() noexcept;

//! Quinoa's own terminate handler
void terminateHandler() noexcept;

// Quinoa's own unexpected handler
void unexpectedHandler () noexcept;

} // namespace quinoa

#endif // Handler_h
