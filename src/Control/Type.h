//******************************************************************************
/*!
  \file      src/Control/Type.h
  \author    J. Bakosi
  \date      Mon Feb  4 16:24:05 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for control and parsing
  \details   Types for control and parsing
*/
//******************************************************************************
#ifndef Type_h
#define Type_h

#include<string>
#include<tuple>

#include<QuinoaTypes.h>

using namespace std;

namespace Quinoa {

namespace control {

//! Physics (methods: collection of models) types
enum PhysicsType { NO_PHYSICS=0,
                   HOMOGENEOUS_DIRICHLET,
                   HOMOGENEOUS_GENERALIZED_DIRICHLET,
                   SPINSFLOW,
                   NUM_PHYSICS
};

//! Hydrodynamics model types
enum HydroType { NO_HYDRO=0,
                 SLM,
                 GLM,
                 NUM_HYDRO
};

//! Material mix model types
enum MixType { NO_MIX=0,
               IEM,
               IECM,
               DIRICHLET,
               GENERALIZED_DIRICHLET,
               NUM_MIX
};

enum BundlePosition { TITLE=0,
                      PHYSICS,
                      HYDRO,
                      MIX,
                      NSTEP,
                      TERM,
                      DT,
                      NSCALAR,
                      NPAR,
                      ECHO
};

// Storage bundle for parsed data
using Bundle = tuple< string,        //!< 0: Title
                      PhysicsType,   //!< 1: Physics
                      HydroType,     //!< 2: Hydrodynamics model
                      MixType,       //!< 3: Material mix model
                      int,           //!< 4: Number of time steps to take
                      real,          //!< 5: Time to terminate time stepping
                      real,          //!< 6: Size of time step
                      int,           //!< 7: Number of mixing scalars
                      int,           //!< 8: Total number of particles
                      int >;         //!< 9: One-liner info every few time steps

} // namespace control

} // namespace Quinoa

#endif // Type_h
