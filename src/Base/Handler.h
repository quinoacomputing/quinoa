//******************************************************************************
/*!
  \file      src/Base/Handler.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:07:41 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Handler functions
  \details   Handler functions
*/
//******************************************************************************
#ifndef Handler_h
#define Handler_h

#include <Exception.h>

namespace tk {

//! Process an exception
ErrCode processException() noexcept;

//! Quinoa's own new handler
void newHandler() noexcept;

//! Quinoa's own terminate handler
void terminateHandler() noexcept;

// Quinoa's own unexpected handler
void unexpectedHandler () noexcept;

} // tk

#endif // Handler_h
