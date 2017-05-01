// *****************************************************************************
/*!
  \file      src/Walker/Integrator.C
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

#include "Integrator.h"
#include "Collector.h"

namespace walker {

extern std::vector< DiffEq > g_diffeqs;

}

using walker::Integrator;

Integrator::Integrator( CProxy_Distributor& hostproxy,
                        CProxy_Collector& collproxy,
                        uint64_t npar ) :
  m_hostproxy( hostproxy ),
  m_collproxy( collproxy ),
  m_particles( npar, g_inputdeck.get< tag::component >().nprop() ),
  m_stat( m_particles,
          g_inputdeck.get< tag::component >().offsetmap( 
            g_inputdeck.depvars() ),
          g_inputdeck.get< tag::stat >(),
          g_inputdeck.get< tag::pdf >(),
          g_inputdeck.get< tag::discr, tag::binsize >() )
// *****************************************************************************
// Constructor
//! \param[in] hostproxy Host proxy to call back to
//! \param[in] collproxy Collector proxy to send results to
//! \param[in] npar Number of particles this integrator advances
//! \author J. Bakosi
// *****************************************************************************
{
  // register with the local branch of the statistics collector
  m_collproxy.ckLocalBranch()->checkin();
  // Tell the Charm++ runtime system to call back to
  // Distributor::registered() once all Integrator chares have registered
  // themselves, i.e., checked in, with their local branch of the statistics
  // merger group, Collector. The reduction is done via creating a callback
  // that invokes the typed reduction client, where m_hostproxy is the proxy
  // on which the reduction target method, registered(), is called upon
  // completion of the reduction.
  contribute(
    CkCallback(CkReductionTarget( Distributor, registered ), m_hostproxy) );
}

void
Integrator::setup( tk::real dt,
                   tk::real t,
                   uint64_t it,
                   const std::map< tk::ctr::Product, tk::real >& moments )
// *****************************************************************************
// Perform setup: set initial conditions and advance a time step
//! \param[in] dt Size of time step
//! \param[in] t Physical time
//! \param[in] it Iteration count
//! \param[in] moments Map of statistical moments
// *****************************************************************************
{
  ic();                           // set initial conditions for all equations
  advance( dt, t, it, moments );  // start time stepping all equations
}

void
Integrator::ic()
// *****************************************************************************
// Set initial conditions
//! \author J. Bakosi
// *****************************************************************************
{
  for (const auto& eq : g_diffeqs) eq.initialize( CkMyPe(), m_particles );
}

void
Integrator::advance( tk::real dt,
                     tk::real t,
                     uint64_t it,
                     const std::map< tk::ctr::Product, tk::real >& moments )
// *****************************************************************************
// Advance all particles owned by this integrator
//! \param[in] dt Size of time step
//! \param[in] t Physical time
//! \param[in] it Iteration count
//! \param[in] moments Map of statistical moments
//! \author J. Bakosi
// *****************************************************************************
{
  // Advance all equations one step in time. At the 0th iteration skip advance
  // but estimate statistics and (potentially) PDFs (at the interval given by
  // the user).
  if (it > 0)
    for (const auto& e : g_diffeqs)
      e.advance( m_particles, CkMyPe(), dt, t, moments );

  if (!g_inputdeck.stat()) {// if no stats to estimate, skip to end of time step
    contribute(
      CkCallback(CkReductionTarget( Distributor, nostat ), m_hostproxy) );
  } else {
    // Accumulate sums for ordinary moments (every time step)
    accumulateOrd();
    // Accumulate sums for ordinary PDFs at select times
    if ( g_inputdeck.pdf() &&
         !((it+1) % g_inputdeck.get< tag::interval, tag::pdf >()) )
      accumulateOrdPDF();
  }
}

void
Integrator::accumulateOrd()
// *****************************************************************************
// Accumulate sums for ordinary moments
//! \author J. Bakosi
// *****************************************************************************
{
  // Accumulate partial sums for ordinary moments
  m_stat.accumulateOrd();
  // Send accumulated ordinary moments to collector for estimation
  m_collproxy.ckLocalBranch()->chareOrd( m_stat.ord() );
}

void
Integrator::accumulateCen( const std::vector< tk::real >& ord )
// *****************************************************************************
// Accumulate sums for central moments
//! \param[in] ord Estimated ordinary moments (collected from all PEs)
//! \author J. Bakosi
// *****************************************************************************
{
  // Accumulate partial sums for central moments
  m_stat.accumulateCen( ord );
  // Send accumulated central moments to host for estimation
  m_collproxy.ckLocalBranch()->chareCen( m_stat.ctr() );
}

void
Integrator::accumulateOrdPDF()
// *****************************************************************************
// Accumulate sums for ordinary PDFs
//! \author J. Bakosi
// *****************************************************************************
{
  // Accumulate partial sums for ordinary PDFs
  m_stat.accumulateOrdPDF();
  // Send accumulated ordinary PDFs to host for estimation
  m_collproxy.ckLocalBranch()->chareOrdPDF( m_stat.oupdf(),
                                            m_stat.obpdf(),
                                            m_stat.otpdf() );
}

void
Integrator::accumulateCenPDF( const std::vector< tk::real >& ord )
// *****************************************************************************
// Accumulate sums for central PDFs
//! \param[in] ord Estimated ordinary moments (collected from all PEs)
//! \author J. Bakosi
// *****************************************************************************
{
  // Accumulate partial sums for central PDFs
  m_stat.accumulateCenPDF( ord );
  // Send accumulated central PDFs to host for estimation
  m_collproxy.ckLocalBranch()->chareCenPDF( m_stat.cupdf(),
                                            m_stat.cbpdf(),
                                            m_stat.ctpdf() );
}

#include "NoWarning/integrator.def.h"
