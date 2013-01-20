//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Sun 20 Jan 2013 07:46:36 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base, select model
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

namespace Quinoa {

//! Available Physics
enum class PhysicsType { HOMOGENEOUS_DIRICHLET,     //!< Homogeneous Dirichlet
                         HOMOGENEOUS_GENDIRICHLET   //!< Hom. Generalized Dir.
};

} // namespace Quinoa

#endif // Control_h
