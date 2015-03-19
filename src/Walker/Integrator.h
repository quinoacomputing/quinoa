//******************************************************************************
/*!
  \file      src/Walker/Integrator.h
  \author    J. Bakosi
  \date      Thu 19 Mar 2015 11:43:50 AM MDT
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

#include <ParticleProperties.h>
#include <DiffEq.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <walker.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <integrator.decl.h>
#include <Walker/InputDeck/InputDeck.h>
#include <Statistics.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::vector< DiffEq > g_diffeqs;

//! Integrator Charm++ chare used to advance differential equations in time
template< class Proxy >
class Integrator : public CBase_Integrator< Proxy > {

  public:
    //! Constructor
    //! \param[in] proxy Host proxy to call back to (here: Distributor)
    //! \param[in] npar Number of particles this integrator advances
    //! \param[in] dt Size of time step
    //! \param[in] it Iteration count
    //! \param[in] moments Map of statistical moments
    explicit Integrator( Proxy& proxy,
                         uint64_t npar,
                         tk::real dt,
                         uint64_t it,
                         const std::map< tk::ctr::Product, tk::real >& moments )
      : m_proxy( proxy ),
        m_particles( npar, g_inputdeck.get< tag::component >().nprop() ),
        m_stat( m_particles,
                g_inputdeck.get< tag::component >().offsetmap( 
                  g_inputdeck.depvars() ),
                g_inputdeck.get< tag::stat >(),
                g_inputdeck.get< tag::pdf >(),
                g_inputdeck.get< tag::discr, tag::binsize >() )
    {
      ic();                 // set initial conditions for all equations
      advance( dt, it, moments );    // start time stepping all equations
    }

    //! Set initial conditions
    void ic() {
      for (const auto& eq : g_diffeqs) eq.initialize( m_particles );
      m_proxy.init();   // signal to host that initialization is complete
    }

    //! Advance all particles owned by this integrator
    //! \param[in] dt Size of time step
    //! \param[in] it Iteration count
    //! \param[in] moments Map of statistical moments
    void advance( tk::real dt,
                  uint64_t it,
                  const std::map< tk::ctr::Product, tk::real >& moments )
    {
      //! Advance all equations one step in time
      for (const auto& e : g_diffeqs)
        e.advance( m_particles, CkMyPe(), dt, moments );
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
    //! \param[in] ord Estimated ordinary moments (collected from all PEs)
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
    //! \param[in] ord Estimated ordinary moments (collected from all PEs)
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

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include <integrator.def.h>
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // Integrator_h
