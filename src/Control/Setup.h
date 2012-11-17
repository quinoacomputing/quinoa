//******************************************************************************
/*!
  \file      src/Control/Setup.h
  \author    J. Bakosi
  \date      Fri 16 Nov 2012 08:49:48 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model setup
  \details   Model setup, see Control.h for available choices
*/
//******************************************************************************
#ifndef Setup_h
#define Setup_h

#include <limits>

#include <Control.h>

namespace Quinoa {

//! Model
const ModelType MODEL_TYPE = ModelType::HOMOGENEOUS_DIRICHLET;

//! Number of scalars
const int NSCALAR = 3;

//! Number of particles
const int NPAR = 10000;

//! Maximum time to simulate
const real TIME = 1.0;

//! Maximum number of time steps to take
const int NSTEP = numeric_limits<int>::max();

} // namespace Quinoa

#endif // Setup_h
