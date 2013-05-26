//******************************************************************************
/*!
  \file      src/Control/BackAssociate.h
  \author    J. Bakosi
  \date      Sun 26 May 2013 05:41:13 PM MDT
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
  const string PhysicsKeyword[] = {
    "no_physics",
    "hommix",
    "homhydro",
    "homrt",
    "spinsflow"
  };
  // PhysicsType -> name
  const string PhysicsName[] = {
    "No physics",
    "Homogeneous material mixing",
    "Homogeneous hydrodynamics",
    "Homogeneous Rayleigh-Taylor",
    "Standalone-Particle Incompressible Navier-Stokes Flow"
  };

  // PositionType -> keyword
  const string PositionKeyword[] = {
    "no_position",
    "invpos",
    "vispos"
  };
  // PositionType -> name
  const string PositionName[] = {
    "No model",
    "Inviscid",
    "Viscous"
  };

  // MassType -> keyword
  const string MassKeyword[] = {
    "no_mass",
    "beta"
  };
  // MassType -> name
  const string MassName[] = {
    "No model",
    "Beta"
  };

  // HydroType -> keyword
  const string HydroKeyword[] = {
    "no_hydro",
    "slm",
    "glm"
  };
  // HydroType -> name
  const string HydroName[] = {
    "No model",
    "Simplified Langevin",
    "Generalized Langevin"
  };

  // MixType -> keyword
  const string MixKeyword[] = {
    "no_mix",
    "iem",
    "iecm",
    "dir",
    "gendir"
  };
  // MixType -> name
  const string MixName[] = {
    "No model",
    "Interaction by exchange with the mean",
    "Interaction by exchange with the conditional mean",
    "Dirichlet",
    "Generalized Dirichlet"
  };

} // namespace associate

} // namespace Quinoa

#endif // BackAssociate_h
