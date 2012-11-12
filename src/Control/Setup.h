//******************************************************************************
/*!
  \file      src/Control/Setup.h
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 10:11:48 AM MST
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
const ModelType MODEL = ModelType::DIRICHLET;

//! Select mix model
const MixModelType MIX_MODEL = MixModelType::DIRICHLET;

//! Select number of scalars
const int NUM_SCALARS = 0;

//! Select velocity model
const VelocityModelType VELOCITY_MODEL = VelocityModelType::NONE;

} // namespace Quinoa

#endif // Setup_h
