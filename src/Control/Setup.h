//******************************************************************************
/*!
  \file      src/Control/Setup.h
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 07:58:18 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model setup
  \details   Model setup, see Control.h for available choices
*/
//******************************************************************************
#ifndef Setup_h
#define Setup_h

#include <Control.h>

namespace Quinoa {

//! Select model
const ModelType MODEL_TYPE = ModelType::HOMOGENEOUS_DIRICHLET;

//! Select number of scalars
const int NSCALAR = 3;

//! Select number of particles per element
const int NPEL = 10;

} // namespace Quinoa

#endif // Setup_h
