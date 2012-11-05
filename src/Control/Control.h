//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Sun 04 Nov 2012 06:52:33 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base, decides which sub-control category is used
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

namespace Quinoa {

//! Mix models
enum MixModel { DIRICHLET=0,            //!< Dirichlet
                GENERALIZED_DIRICHLET,  //!< Generalized Dirichlet
                NUM_MIX_MODELS
};

//! Velocity models
enum VelocityModel { SIMPLIFIED_LANGEVIN=0,//!< Homogeneous Simplified Langevin
                     GENERALIZED_LANGEVIN, //!< Homogeneous Generalized Langevin
                     NUM_VELOCITY_MODELS
};

} // namespace Quinoa

#endif // Control_h
