//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 10:59:06 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base, select model
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

namespace Quinoa {

//! Available physics (methods) options
enum class PhysicsType {
  HOMOGENEOUS_DIRICHLET,     //!< Homogeneous Dirichlet
  HOMOGENEOUS_GENDIRICHLET,  //!< Homogeneous Generalized Dirichlet
  SPINSFLOW           //!< Standalone-Particle Incompressible Navier-Stokes Flow
};

//! Available hydrodynamics models
enum class HydroType {
  SLM,     //!< Simplified Langevin model
  GLM      //!< Generalized Langevin model
};

} // namespace Quinoa

#endif // Control_h
