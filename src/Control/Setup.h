//******************************************************************************
/*!
  \file      src/Control/Setup.h
  \author    J. Bakosi
  \date      Sun 04 Nov 2012 07:14:03 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Basic setup: specify models, parameters, etc.
  \details   Basic setup: specify models, parameters, etc.
*/
//******************************************************************************
#ifndef Setup_h
#define Setup_h

#include <Control.h>

namespace Quinoa {

//! Select mix model; see Control.h for available mix models
const MixModel mixModel = DIRICHLET;

} // namespace Quinoa

#endif // Setup_h
