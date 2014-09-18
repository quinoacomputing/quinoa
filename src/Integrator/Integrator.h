//******************************************************************************
/*!
  \file      src/Integrator/Integrator.h
  \author    J. Bakosi
  \date      Tue 16 Sep 2014 07:29:47 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Integrator used to advance ordinary and stochastic differential eqs.
  \details   Integrator used to advance ordinary and stochastic differential
             equations.
*/
//******************************************************************************
#ifndef Integrator_h
#define Integrator_h

#include <ParticleProperties.h>
#include <DiffEq.h>
#include <quinoa.decl.h>
#include <integrator.decl.h>
#include <Quinoa/InputDeck/InputDeck.h>
#include <Statistics.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck;
extern std::vector< DiffEq > g_diffeqs;

//! Integrator used to advance ODEs and SDEs in time
template< class Proxy >
class Integrator : public CBase_Integrator< Proxy > {

  public:
    //! Constructor
    explicit Integrator( Proxy& proxy, uint64_t npar, tk::real dt ) :
      m_proxy( proxy ),
      m_particles( npar, g_inputdeck.get< tag::component >().nprop() ),
      m_stat( m_particles )
    {
      ic();             // set initial conditions
      advance( dt );    // start time stepping
    }

    //! Set initial conditions
    void ic() {
      for (const auto& eq : g_diffeqs) eq.initialize( m_particles );
      m_proxy.init();   // signal to host that initialization is complete
    }

    //! Advance all equations one step in time
    void advance( tk::real dt ) {
      for (const auto& e : g_diffeqs) e.advance( m_particles, CkMyPe(), dt );
      accumulateOrd();    // start accumulating ordinary moments
    }

    // Accumulate sums for central moments
    void accumulateCen( const std::vector< tk::real >& ord ) {
      // Accumulate partial sums for central moments
      m_stat.accumulateCen( ord );
      // Send accumulated central moments to host for estimation
      m_proxy.estimateCen( m_stat.ctr() );
    }

    // Accumulate sums for PDFs
    void accumulatePDF() {
      // Accumulate partial sums for univariate PDFs
      m_stat.accumulatePDF();
      // Send accumulated PDFs to host for estimation
      m_proxy.estimatePDF( m_stat.updf(), m_stat.bpdf() );
    }

  private:
    // Accumulate sums for ordinary moments
    void accumulateOrd() {
      // Accumulate partial sums for ordinary moments
      m_stat.accumulateOrd();
      // Send accumulated ordinary moments to host for estimation
      m_proxy.estimateOrd( m_stat.ord() );
    }

    Proxy m_proxy;              //!< Host proxy
    ParProps m_particles;       //!< Particle properties operated on
    Statistics m_stat;          //!< Statistics
};

} // quinoa::

#define CK_TEMPLATES_ONLY
#include <integrator.def.h>
#undef CK_TEMPLATES_ONLY

#endif // Integrator_h
