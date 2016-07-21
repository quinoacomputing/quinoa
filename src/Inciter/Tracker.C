// *****************************************************************************
/*!
  \file      src/Inciter/Tracker.C
  \author    F.J. Gonzalez
  \date      Thu 21 Jul 2016 02:48:27 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Tracker advances Lagrangian particles passively tracking fluid flow
  \details   Tracker advances Lagrangian particles passively tracking fluid
    flow. There are a potentially large number of Tracker Charm++ chares
    created by Conductor and Partitioner. Each tracker gets a chunk of the full
    load (part of the total number particles) and does the same: initializes and
    advances a particles in time and space (across a mesh) based on a velocity
    field.
*/
// *****************************************************************************

#include "Tracker.h"
#include "NoWarning/performer.decl.h"

#include "tracker.def.h"
