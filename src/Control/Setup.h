//******************************************************************************
/*!
  \file      src/Control/Setup.h
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 08:39:41 PM MST
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
const MixModel MIX_MODEL = MIX_DIRICHLET;

//! Select velocity model; see Control.h for available velocity models
const VelocityModel VELOCITY_MODEL = VEL_NONE;

} // namespace Quinoa

#endif // Setup_h
