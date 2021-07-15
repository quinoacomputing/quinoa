// *****************************************************************************
/*!
  \file      src/Walker/Integrator.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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

#include "Integrator.hpp"
#include "Collector.hpp"

namespace walker {

extern std::vector< DiffEq > g_diffeqs;

}

using walker::Integrator;

Integrator::Integrator( CProxy_Distributor hostproxy,
                        CProxy_Collector collproxy,
                        tk::CProxy_ParticleWriter particlewriterproxy,
                        uint64_t npar ) :
  m_host( hostproxy ),
  m_coll( collproxy ),
  m_particlewriter( particlewriterproxy ),
  m_particles( npar, g_inputdeck.get< tag::component >().nprop() ),
  m_stat( m_particles,
          g_inputdeck.get< tag::component >().offsetmap( g_inputdeck ),
          g_inputdeck.get< tag::stat >(),
          g_inputdeck.get< tag::pdf >(),
          g_inputdeck.get< tag::discr, tag::binsize >() ),
  m_dt( 0.0 ),
  m_t( 0.0 ),
  m_it( 0 ),
  m_itp( 0 )
// *****************************************************************************
// Constructor
//! \param[in] hostproxy Host proxy to call back to
//! \param[in] collproxy Collector proxy to send results to
//! \param[in] particlewriterproxy Particle writer proxy to use
//! \param[in] npar Number of particles this integrator advances
// *****************************************************************************
{
  // register with the local branch of the statistics collector
  m_coll.ckLocalBranch()->checkin();
  // Tell the Charm++ runtime system to call back to
  // Distributor::registered() once all Integrator chares have registered
  // themselves, i.e., checked in, with their local branch of the statistics
  // merger group, Collector. The reduction is done via creating a callback
  // that invokes the typed reduction client, where m_host is the proxy
  // on which the reduction target method, registered(), is called upon
  // completion of the reduction.
  contribute( CkCallback(CkReductionTarget(Distributor, registered), m_host) );
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
// *****************************************************************************
{
  // Advance all equations one step in time. At the 0th iteration skip advance
  // but estimate statistics and (potentially) PDFs (at the interval given by
  // the user).
  if (it > 0)
    for (const auto& e : g_diffeqs)
      e.advance( m_particles, CkMyPe(), dt, t, moments );

  // Save time stepping data
  m_dt = dt;
  m_t = t;
  m_it = it;

  // Contribute number of particles we hit the particles output frequency
  auto poseq =
    !g_inputdeck.get< tag::param, tag::position, tag::depvar >().empty();
  const auto parfreq = g_inputdeck.get< tag::interval, tag::particles >();

  CkCallback c( CkIndex_Integrator::out(), thisProxy[thisIndex] );

  if (poseq && !((m_it+1) % parfreq))
    m_particlewriter[ CkMyNode() ].npar( m_particles.nunk(), c );
  else
    c.send();
}

void
Integrator::out()
// *****************************************************************************
// Output particle positions to file
// *****************************************************************************
{
  auto poseq =
    !g_inputdeck.get< tag::param, tag::position, tag::depvar >().empty();
  const auto parfreq = g_inputdeck.get< tag::interval, tag::particles >();

  CkCallback c( CkIndex_Integrator::accumulate(), thisProxy[thisIndex] );

  // Output particles data to file if we hit the particles output frequency
  if (poseq && !((m_it+1) % parfreq)) {
    // query position eq offset in particle array (0: only first particle pos)
    auto po = g_inputdeck.get< tag::component >().offset< tag::position >( 0 );
    // output particle positions to file
    m_particlewriter[ CkMyNode() ].
      writeCoords( m_itp++, m_particles.extract(0,po),
        m_particles.extract(1,po), m_particles.extract(2,po), c );
  } else {
    c.send();
  }
}

void
Integrator::accumulate()
// *****************************************************************************
// Start collecting statistics
// *****************************************************************************
{
  if (!g_inputdeck.stat()) {// if no stats to estimate, skip to end of time step
    contribute( CkCallback(CkReductionTarget(Distributor, nostat), m_host) );
  } else {
    // Accumulate sums for ordinary moments (every time step)
    accumulateOrd( m_it, m_t, m_dt );
  }
}

void
Integrator::accumulateOrd( uint64_t it, tk::real t, tk::real dt )
// *****************************************************************************
// Accumulate sums for ordinary moments and ordinary PDFs
//! \param[in] it Iteration count
//! \param[in] t Physical time
//! \param[in] dt Time step size
// *****************************************************************************
{
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto pdffreq = g_inputdeck.get< tag::interval, tag::pdf >();

  // Accumulate partial sums for ordinary moments
  m_stat.accumulateOrd();
  // Accumulate sums for ordinary PDFs at first and last iterations and at
  // select times
  if ( g_inputdeck.pdf() &&
       ( it == 0 ||
         !((it+1) % pdffreq) ||
         (std::fabs(t+dt-term) < eps && (it+1) >= nstep) ) )
    m_stat.accumulateOrdPDF();

  // Send accumulated ordinary moments and ordinary PDFs to collector for
  // estimation
  m_coll.ckLocalBranch()->chareOrd( m_stat.ord(),
                                    m_stat.oupdf(),
                                    m_stat.obpdf(),
                                    m_stat.otpdf() );
}

void
Integrator::accumulateCen( uint64_t it,
                           tk::real t,
                           tk::real dt,
                           const std::vector< tk::real >& ord )
// *****************************************************************************
// Accumulate sums for central moments and central PDFs
//! \param[in] it Iteration count
//! \param[in] t Physical time
//! \param[in] dt Time step size
//! \param[in] ord Estimated ordinary moments (collected from all PEs)
// *****************************************************************************
{
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto pdffreq = g_inputdeck.get< tag::interval, tag::pdf >();

  // Accumulate partial sums for central moments
  m_stat.accumulateCen( ord );
  // Accumulate partial sums for central PDFs at first and last iteraions and
  // at select times
  if ( g_inputdeck.pdf() &&
       ( it == 0 ||
         !((it+1) % pdffreq) ||
         (std::fabs(t+dt-term) < eps && (it+1) >= nstep) ) )
    m_stat.accumulateCenPDF( ord );

  // Send accumulated central moments to host for estimation
  m_coll.ckLocalBranch()->chareCen( m_stat.ctr(),
                                    m_stat.cupdf(),
                                    m_stat.cbpdf(),
                                    m_stat.ctpdf() );
}

#include "NoWarning/integrator.def.h"
