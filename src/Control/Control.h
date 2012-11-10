//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Sat 10 Nov 2012 07:56:35 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base, select model
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

namespace Quinoa {

//! Models
enum class ModelType { DIRICHLET,              //!< Dirichlet
                       GENERALIZED_DIRICHLET   //!< Generalized Dirichlet
};

//! Mix models
enum class MixModelType { NONE,                   //!< No mix model
                          DIRICHLET,              //!< Dirichlet
                          GENERALIZED_DIRICHLET   //!< Generalized Dirichlet
};

//! Velocity models
enum class VelocityModelType { NONE,                 //!< No velocity model
                               SIMPLIFIED_LANGEVIN,  //!< Simplified Langevin
                               GENERALIZED_LANGEVIN  //!< Generalized Langevin
};

} // namespace Quinoa

#endif // Control_h
