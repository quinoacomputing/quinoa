// *****************************************************************************
/*!
  \file      src/Walker/Integrator.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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

#include "Types.h"
#include "Tags.h"
#include "StatCtr.h"
#include "DiffEq.h"
#include "Particles.h"
#include "SystemComponents.h"
#include "Statistics.h"
#include "Walker/InputDeck/InputDeck.h"

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
    explicit Integrator( CProxy_Distributor& hostproxy,
                         CProxy_Collector& collproxy,
                         uint64_t npar );

    //! Migrate constructor
    explicit Integrator( CkMigrateMessage* ) :
      // WARNING: This is a "blind" copy of the standard constructor initializer
      // list - it must be changed for migration to be correct.      
      m_hostproxy(),
      m_collproxy(),
      m_particles( 0, g_inputdeck.get< tag::component >().nprop() ),
      m_stat( m_particles,
                g_inputdeck.get< tag::component >().offsetmap( 
                  g_inputdeck.depvars() ),
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

    // Accumulate sums for ordinary moments
    void accumulateOrd();

    // Accumulate sums for central moments
    void accumulateCen( const std::vector< tk::real >& ord );

    // Accumulate sums for ordinary PDFs
    void accumulateOrdPDF();

    // Accumulate sums for central PDFs
    void accumulateCenPDF( const std::vector< tk::real >& ord );

  private:
    CProxy_Distributor m_hostproxy;     //!< Host proxy
    CProxy_Collector m_collproxy;       //!< Collector proxy
    tk::Particles m_particles;          //!< Particle properties
    tk::Statistics m_stat;              //!< Statistics
};

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

} // walker::

#endif // Integrator_h
