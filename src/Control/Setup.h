//******************************************************************************
/*!
  \file      src/Control/Setup.h
  \author    J. Bakosi
  \date      Fri Nov 16 08:49:49 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model setup
  \details   Model setup, see Control.h for available choices
*/
//******************************************************************************
#ifndef Setup_h
#define Setup_h

#include <Control.h>

namespace Quinoa {

//! Model
const ModelType MODEL_TYPE = ModelType::HOMOGENEOUS_DIRICHLET;

//! Number of scalars
const int NSCALAR = 3;

//! Number of particles
const int NPAR = 10;

} // namespace Quinoa

#endif // Setup_h
