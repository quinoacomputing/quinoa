//******************************************************************************
/*!
  \file      src/Walker/Integrator.h
  \author    J. Bakosi
  \date      Thu 22 Jan 2015 06:12:56 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Integrator used to advance ordinary and stochastic differential eqs.
  \details   Integrator used to advance ordinary and stochastic differential
             equations.
*/
//******************************************************************************
#ifndef Integrator_h
#define Integrator_h

#include <ParticleProperties.h>
#include <DiffEq.h>
#include <walker.decl.h>
#include <integrator.decl.h>
#include <Walker/InputDeck/InputDeck.h>
#include <Statistics.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::vector< DiffEq > g_diffeqs;

//! Integrator used to advance ODEs and SDEs in time
template< class Proxy >
class Integrator : public CBase_Integrator< Proxy > {

  public:
    //! Constructor
    explicit Integrator( Proxy& proxy, uint64_t npar, tk::real dt, uint64_t it )
      : m_proxy( proxy ),
        m_particles( npar, g_inputdeck.get< tag::component >().nprop() ),
        m_stat( m_particles,
                g_inputdeck.get< tag::component >().offsetmap( 
                  g_inputdeck.depvars() ),
                g_inputdeck.get< tag::stat >(),
                g_inputdeck.get< tag::pdf >(),
                g_inputdeck.get< tag::discr, tag::binsize >() )
    {
      ic();                 // set initial conditions
      advance( dt, it );    // start time stepping
    }

    //! Set initial conditions
    void ic() {
      for (const auto& eq : g_diffeqs) eq.initialize( m_particles, m_stat );
      m_proxy.init();   // signal to host that initialization is complete
    }

    //! Advance all particles owned by this integrator
    void advance( tk::real dt, uint64_t it ) {
      //! Advance all equations one step in time
      for (const auto& e : g_diffeqs) e.advance( m_particles, CkMyPe(), dt );
      // Accumulate sums for ordinary moments (every time step)
      accumulateOrd();
      // Accumulate sums for ordinary PDFs at select times
      if ( !(it % g_inputdeck.get< tag::interval, tag::pdf >()) )
        accumulateOrdPDF();
    }

    // Accumulate sums for ordinary moments
    void accumulateOrd() {
      // Accumulate partial sums for ordinary moments
      m_stat.accumulateOrd();
      // Send accumulated ordinary moments to host for estimation
      m_proxy.estimateOrd( m_stat.ord() );
    }

    // Accumulate sums for central moments
    void accumulateCen( const std::vector< tk::real >& ord ) {
      // Accumulate partial sums for central moments
      m_stat.accumulateCen( ord );
      // Send accumulated central moments to host for estimation
      m_proxy.estimateCen( m_stat.ctr() );
    }

    // Accumulate sums for ordinary PDFs
    void accumulateOrdPDF() {
      // Accumulate partial sums for ordinary PDFs
      m_stat.accumulateOrdPDF();
      // Send accumulated ordinary PDFs to host for estimation
      m_proxy.estimateOrdPDF( m_stat.oupdf(), m_stat.obpdf(), m_stat.otpdf() );
    }

    // Accumulate sums for central PDFs
    void accumulateCenPDF( const std::vector< tk::real >& ord ) {
      // Accumulate partial sums for central PDFs
      m_stat.accumulateCenPDF( ord );
      // Send accumulated central PDFs to host for estimation
      m_proxy.estimateCenPDF( m_stat.cupdf(), m_stat.cbpdf(), m_stat.ctpdf() );
    }

  private:
    Proxy m_proxy;                      //!< Host proxy
    tk::ParProps m_particles;           //!< Particle properties
    tk::Statistics m_stat;              //!< Statistics
};

} // walker::

#define CK_TEMPLATES_ONLY
#include <integrator.def.h>
#undef CK_TEMPLATES_ONLY

#endif // Integrator_h
