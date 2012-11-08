//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Thu 08 Nov 2012 05:40:47 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base, select model
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

namespace Quinoa {

//! Models
enum ModelType { DIRICHLET=0,            //!< Dirichlet
                 GENERALIZED_DIRICHLET,  //!< Generalized Dirichlet
                 NUM_MODELS
};

//! Mix models
enum MixModelType { MIX_NONE=0,                 //!< No mix model
                    MIX_DIRICHLET,              //!< Dirichlet
                    MIX_GENERALIZED_DIRICHLET,  //!< Generalized Dirichlet
                    NUM_MIX_MODELS
};

//! Velocity models
enum VelocityModelType { VEL_NONE=0,                //!< No velocity model
                         VEL_SIMPLIFIED_LANGEVIN,   //!< Simplified Langevin
                         VEL_GENERALIZED_LANGEVIN,  //!< Generalized Langevin
                         NUM_VEL_MODELS
};

} // namespace Quinoa

#endif // Control_h
