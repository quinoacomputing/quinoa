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
  m_conductor( conductor ),
  m_performer( performer )
// *****************************************************************************
//  Constructor
//! \param[in] conductor Host (Conductor) proxy
//! \author F.J. Gonzalez
// *****************************************************************************
{
  /*auto& x = m_pcoord[0];
  auto& y = m_pcoord[1];
  auto& z = m_pcoord[2];
  for (std::size_t i=0; i<10; ++i) {
    tk::real a=i;
    x.push_back(a/10);
    y.push_back(a/10);
    z.push_back(a/10);
  }*/
  m_performer[thisIndex].genPar(10);
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
