// *****************************************************************************
/*!
  \file      src/Inciter/Tracker.h
  \author    F.J. Gonzalez
  \date      Thu 21 Jul 2016 02:53:44 PM MDT
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
#ifndef Tracker_h
#define Tracker_h

#include <array>
#include <iostream>
#include <unordered_map>

#include "Types.h"
#include "Macro.h"

#include "NoWarning/conductor.decl.h"
#include "NoWarning/tracker.decl.h"

namespace inciter {

//! \brief Tracker Charm++ chare used to advance Lagrangian particles tracking
//!   fluid flow in space and time
template< class WorkerProxy >
class Tracker : public CBase_Tracker< WorkerProxy > {

  private:
    using ConductorProxy = CProxy_Conductor;
    using Array = CBase_Tracker< WorkerProxy >;

  public:
    //! Constructor
    explicit Tracker( const ConductorProxy& conductor,
                      const WorkerProxy& worker ) :
      m_conductor( conductor ),
      m_worker( worker ),
      m_pcoord() {}

    //! Migrate constructor
    explicit Tracker( CkMigrateMessage* ) :
      m_conductor(), m_worker(), m_pcoord() {}

    //! Advance particles to next stage
    void advance( tk::real dt, uint64_t it, tk::real t ) {
      std::cout << "tracker id " << Array::thisIndex << ": advance\n";
      m_worker[ Array::thisIndex ].recPartLoc( m_pcoord[0],
                                               m_pcoord[1],
                                               m_pcoord[2] );
      IGNORE(dt);
      IGNORE(it);
      IGNORE(t);
    }

  private:
    ConductorProxy m_conductor;                          //!< Conductor proxy
    WorkerProxy m_worker;                                //!< Worker proxy
    std::array< std::vector< tk::real >, 3 > m_pcoord;   //!< Particle coords
};

} // inciter::

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wreorder"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wreorder"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#define CK_TEMPLATES_ONLY
#include "tracker.def.h"
#undef CK_TEMPLATES_ONLY

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // Tracker_h
