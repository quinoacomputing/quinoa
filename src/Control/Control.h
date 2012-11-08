//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 08:38:53 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base, decides which sub-control category is used
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

namespace Quinoa {

//! Mix models
enum MixModel { MIX_NONE=0,                 //!< No mix model
                MIX_DIRICHLET,              //!< Dirichlet
                MIX_GENERALIZED_DIRICHLET,  //!< Generalized Dirichlet
                NUM_MIX_MODELS
};

//! Velocity models
enum VelocityModel { VEL_NONE=0,                //!< No velocity model
                     VEL_SIMPLIFIED_LANGEVIN,   //!< Simplified Langevin
                     VEL_GENERALIZED_LANGEVIN,  //!< Generalized Langevin
                     NUM_VEL_MODELS
};

} // namespace Quinoa

#endif // Control_h
