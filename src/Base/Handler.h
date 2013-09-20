//******************************************************************************
/*!
  \file      src/Base/Handler.h
  \author    J. Bakosi
  \date      Thu 19 Sep 2013 08:56:25 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Handler functions
  \details   Handler functions
*/
//******************************************************************************
#ifndef Handler_h
#define Handler_h

#include <Exception.h>

namespace quinoa {

//! Process an exception
ErrCode processException() noexcept;

//! Quinoa's own new handler
void newHandler() noexcept;

//! Quinoa's own terminate handler
void terminateHandler() noexcept;

// Quinoa's own unexpected handler
void unexpectedHandler () noexcept;

} // namespace quinoa

#endif // Handler_h
