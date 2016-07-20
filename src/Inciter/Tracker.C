// *****************************************************************************
/*!
  \file      src/Inciter/Tracker.C
  \author    F.J. Gonzalez
  \date      Wed Jul 20 14:27:33 2016
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \details   Tracker advances Lagrangian particles passively tracking fluid
    flow.. There are a potentially large number of Tracker Charm++ chares
    created by Conductor. Each tracker gets a chunk of the full load (part of
    the total number particles) and does the same: initializes and advances a
    particles in time based on a velocity field.
*/
// *****************************************************************************

#include "Tracker.h"
#include "Macro.h"

using inciter::Tracker;

Tracker::Tracker( const ConductorProxy& conductor,
                  const PerformerProxy& performer ) :
  m_conductor( conductor )
// *****************************************************************************
//  Constructor
//! \param[in] conductor Host (Conductor) proxy
//! \author F.J. Gonzalez
// *****************************************************************************
{
}

void
Tracker::advance( tk::real dt, uint64_t it, tk::real t )
// *****************************************************************************
// Advance equations to next stage
//! \param[in] dt Size of time step
//! \param[in] it Iteration count
//! \param[in] t Physical time
//! \author F.J. Gonzalez
// *****************************************************************************
{
  IGNORE(dt);
  IGNORE(it);
  IGNORE(t);
}

#include "NoWarning/tracker.def.h"
