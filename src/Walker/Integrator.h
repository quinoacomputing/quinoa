//******************************************************************************
/*!
  \file      src/Walker/Integrator.h
  \author    J. Bakosi
  \date      Wed 15 Jul 2015 08:56:16 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Integrator advances differential equations
  \details   Integrator advances differential equations. There are a potentially
    large number of Integrator Charm++ chares created by Distributor. Each
    integrator gets a chunk of the full load and does the same: initializes and
    advances multiple ordinary or stochastic differential equations in time.
    Note that there is no spatial dependence, these equations describe spatially
    homogeneous processes.
*/
//******************************************************************************
#ifndef Integrator_h
#define Integrator_h

#include <vector>
#include <map>
#include <cstdint>

#include "Types.h"
#include "Tags.h"
#include "StatCtr.h"
#include "DiffEq.h"
#include "ParticleProperties.h"
#include "SystemComponents.h"
#include "Statistics.h"
#include "Walker/InputDeck/InputDeck.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "integrator.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::vector< DiffEq > g_diffeqs;

//! Integrator Charm++ chare used to advance differential equations in time
class Integrator : public CBase_Integrator {

  public:
    //! Constructor
    explicit Integrator( CProxy_Distributor& hostproxy,
                         CProxy_Collector& collproxy,
                         uint64_t npar );

    //! Migrate constructor
    explicit Integrator( CkMigrateMessage* ) :
      m_stat( m_particles,
                g_inputdeck.get< tag::component >().offsetmap( 
                  g_inputdeck.depvars() ),
                g_inputdeck.get< tag::stat >(),
                g_inputdeck.get< tag::pdf >(),
                g_inputdeck.get< tag::discr, tag::binsize >() ) {}

    //! Perform setup: set initial conditions and advance a time step
    void setup( tk::real dt,
                 uint64_t it,
                 const std::map< tk::ctr::Product, tk::real >& moments );

    //! Set initial conditions
    void ic();

    //! Advance all particles owned by this integrator
    void advance( tk::real dt,
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
    tk::ParProps m_particles;           //!< Particle properties
    tk::Statistics m_stat;              //!< Statistics
    bool m_nostat;                      //!< Any statistics to estimate?
};

} // walker::

#endif // Integrator_h
