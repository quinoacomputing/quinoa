//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 12:17:38 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base, select model
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

namespace Quinoa {

//! Models
enum class ModelType { HOMOGENEOUS_DIRICHLET,     //!< Homogeneous Dirichlet
                       HOMOGENEOUS_GENDIRICHLET   //!< Hom. Generalized Dir.
};

} // namespace Quinoa

#endif // Control_h
