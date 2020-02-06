// *****************************************************************************
/*!
  \file      src/Walker/Integrator.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Integrator advances differential equations
  \details   Integrator advances differential equations. There are a potentially
    large number of Integrator Charm++ chares created by Distributor. Each
    integrator gets a chunk of the full load and does the same: initializes and
    advances multiple ordinary or stochastic differential equations in time.
    Note that there is no spatial dependence, these equations describe spatially
    homogeneous processes.
*/
// *****************************************************************************
#ifndef Integrator_h
#define Integrator_h

#include <vector>
#include <map>
#include <cstdint>

#include "Types.hpp"
#include "Tags.hpp"
#include "StatCtr.hpp"
#include "DiffEq.hpp"
#include "Particles.hpp"
#include "SystemComponents.hpp"
#include "Statistics.hpp"
#include "Walker/InputDeck/InputDeck.hpp"

#include "NoWarning/integrator.decl.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

//! Integrator Charm++ chare used to advance differential equations in time
class Integrator : public CBase_Integrator {

  public:
    //! Constructor
    explicit Integrator( CProxy_Distributor hostproxy,
                         CProxy_Collector collproxy,
                         tk::CProxy_ParticleWriter particlewriterproxy,
                         uint64_t npar );

    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit Integrator( CkMigrateMessage* ) :
      m_particles( 0, g_inputdeck.get< tag::component >().nprop() ),
      m_stat( m_particles,
              g_inputdeck.get< tag::component >().offsetmap( g_inputdeck ),
              g_inputdeck.get< tag::stat >(),
              g_inputdeck.get< tag::pdf >(),
              g_inputdeck.get< tag::discr, tag::binsize >() ) {}

    //! Perform setup: set initial conditions and advance a time step
    void setup( tk::real dt,
                tk::real t,
                uint64_t it,
                const std::map< tk::ctr::Product, tk::real >& moments );

    //! Set initial conditions
    void ic();

    //! Advance all particles owned by this integrator
    void advance( tk::real dt,
                  tk::real t,
                  uint64_t it,
                  const std::map< tk::ctr::Product, tk::real >& moments );

    // Accumulate sums for ordinary moments and ordinary PDFs
    void accumulateOrd( uint64_t it, tk::real t, tk::real dt );

    // Accumulate sums for central moments and central PDFs
    void accumulateCen( uint64_t it,
                        tk::real t,
                        tk::real dt,
                        const std::vector< tk::real >& ord );

  private:
    CProxy_Distributor m_host;     //!< Host proxy
    CProxy_Collector m_coll;       //!< Collector proxy
    tk::CProxy_ParticleWriter m_particlewriter;  //!< Particle writer proxy
    tk::Particles m_particles;     //!< Particle properties
    tk::Statistics m_stat;         //!< Statistics
    uint64_t m_itp;                //!< Particle position output iteration count
};

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

} // walker::

#endif // Integrator_h
