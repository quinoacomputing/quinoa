//******************************************************************************
/*!
  \file      src/Control/BackAssociate.h
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 12:04:43 PM MDT
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

  // PhysicsType -> keyword
  const string PhysicsKeyword[control::NUM_PHYSICS] = {
    "no_physics",
    "hommix",
    "homhydro",
    "spinsflow"
  };
  // PhysicsType -> name
  const string PhysicsName[control::NUM_PHYSICS] = {
    "No physics",
    "Homogeneous material mixing",
    "Homogeneous hydrodynamics",
    "Standalone-Particle Incompressible Navier-Stokes Flow"
  };

  // HydroType -> keyword
  const string HydroKeyword[control::NUM_HYDRO] = {
    "no_hydro",
    "slm",
    "glm"
  };
  // HydroType -> name
  const string HydroName[control::NUM_HYDRO] = {
    "No model",
    "Simplified Langevin",
    "Generalized Langevin"
  };

  // MixType -> keyword
  const string MixKeyword[control::NUM_MIX] = {
    "no_mix",
    "iem",
    "iecm",
    "dir",
    "gendir"
  };
  // MixType -> name
  const string MixName[control::NUM_MIX] = {
    "No model",
    "Interaction by exchange with the mean",
    "Interaction by exchange with the conditional mean",
    "Dirichlet",
    "Generalized Dirichlet"
  };

} // namespace associate

} // namespace Quinoa

#endif // BackAssociate_h
