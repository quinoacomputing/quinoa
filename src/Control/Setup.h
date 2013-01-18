//******************************************************************************
/*!
  \file      src/Control/Setup.h
  \author    J. Bakosi
  \date      Thu 17 Jan 2013 11:26:24 PM MST
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

//! Number of prognostic scalars
const int NSCALAR = 2;

//! Number of particles
const int NPAR = 400000;

//! Maximum time to simulate
const real TIME = 1.0;

//! Maximum number of time steps to take
const int NSTEP = numeric_limits<int>::max();

//! One-liner info in every few time steps
const int ECHO = 10;

} // namespace Quinoa

#endif // Setup_h
