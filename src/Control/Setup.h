//******************************************************************************
/*!
  \file      src/Control/Setup.h
  \author    J. Bakosi
  \date      Sat 02 Feb 2013 08:12:39 AM MST
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

//! Physics
const control::PhysicsType PHYSICS_TYPE = control::PhysicsType::SPINSFLOW;

//! Hydrodynamics model
const control::HydroType HYDRO_TYPE = control::HydroType::SLM;

//! Mesh filename
const string MESH_FILENAME = "cylinder.msh";

//! Number of prognostic scalars
const int NSCALAR = 2;

//! Number of particles
const int NPAR = 4000;

//! Maximum time to simulate
const real TIME = 1.0;

//! Maximum number of time steps to take
const int NSTEP = numeric_limits<int>::max();

//! One-liner info in every few time steps
const int ECHO = 1;

} // namespace Quinoa

#endif // Setup_h
