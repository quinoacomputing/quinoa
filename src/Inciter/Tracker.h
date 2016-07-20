// *****************************************************************************
/*!
  \file      src/Inciter/Tracker.h
  \author    F.J. Gonzalez
  \date      Wed Jul 20 14:26:23 2016
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Tracker advances Lagrangian particles passively tracking fluid flow
  \details   Tracker advances Lagrangian particles passively tracking fluid
    flow.. There are a potentially large number of Tracker Charm++ chares
    created by Conductor. Each tracker gets a chunk of the full load (part of
    the total number particles) and does the same: initializes and advances a
    particles in time based on a velocity field.
*/
// *****************************************************************************
#ifndef Tracker_h
#define Tracker_h

#include <unordered_map>

#include "Types.h"

#include "NoWarning/conductor.decl.h"
#include "NoWarning/tracker.decl.h"

namespace inciter {

//! \brief Tracker Charm++ chare used to advance Lagrangian particles tracking
//!   fluid flow in space and time
class Tracker : public CBase_Tracker {

  private:
    using ConductorProxy = CProxy_Conductor;
    using PerformerProxy = CProxy_Performer;

  public:
    //! Constructor
    explicit Tracker( const ConductorProxy& conductor,
                      const PerformerProxy& performer );

    //! Migrate constructor
    explicit Tracker( CkMigrateMessage* ) : m_conductor() {}

    //! Advance particles to next stage
    void advance( tk::real dt, uint64_t it, tk::real t );

  private:
    ConductorProxy m_conductor;         //!< Conductor proxy
};

} // inciter::

#endif // Tracker_h
