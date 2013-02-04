//******************************************************************************
/*!
  \file      src/Control/BackAssociate.h
  \author    J. Bakosi
  \date      Mon Feb  4 16:29:01 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Backward, data to keyword/name associations
  \details   Backward, data to keyword/name associations
*/
//******************************************************************************
#ifndef BackAssociate_h
#define BackAssociate_h

namespace Quinoa {

namespace associate {

  // Editing anything below should be accompanied by the corresponding changes
  // in FwdAssociate.h as well.

  using namespace control;

  // PhysicsType -> keyword
  const string PhysicsKeyword[NUM_PHYSICS] = {
    "no_physics",
    "homdir",
    "homgendir",
    "spinsflow"
  };
  // PhysicsType -> name
  const string PhysicsName[NUM_PHYSICS] = {
    "No physics selected",
    "Homogeneous Dirichlet",
    "Homogeneous generalized Dirichlet",
    "Standalone-Particle Incompressible Navier-Stokes Flow"
  };

  // HydroType -> keyword
  const string HydroKeyword[NUM_HYDRO] = {
    "no_hydro",
    "slm",
    "glm"
  };
  // HydroType -> name
  const string HydroName[NUM_HYDRO] = {
    "No hydrodynamics model selected",
    "Simplified Langevin model",
    "Generalized Langevin model"
  };

  // MixType -> keyword
  const string MixKeyword[NUM_MIX] = {
    "no_mix",
    "iem",
    "iecm",
    "dir",
    "gendir"
  };
  // MixType -> name
  const string MixName[NUM_MIX] = {
    "No material mix model",
    "Interaction by exchange with the mean",
    "Interaction by exchange with the conditional mean",
    "Dirichlet",
    "Generalized Dirichlet"
  };

} // namespace associate

} // namespace Quinoa

#endif // BackAssociate_h
