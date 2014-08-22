//******************************************************************************
/*!
  \file      src/Integrator/Integrator.h
  \author    J. Bakosi
  \date      Fri 15 Aug 2014 11:50:35 AM MDT
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

namespace quinoa {

extern ctr::InputDeck g_inputdeck;
extern std::vector< DiffEq > g_diffeqs;

//! Integrator used to advance ODEs and SDEs in time
template< class Proxy >
class Integrator : public CBase_Integrator< Proxy > {

  public:
    //! Constructor
    explicit Integrator( Proxy& proxy, uint64_t npar ) :
      m_proxy( proxy ),
      m_particles( npar, g_inputdeck.get< tag::component >().nprop() ),
      m_it( 0 ),
      m_t( 0.0 ),
      m_dt( g_inputdeck.get< tag::discr, tag::dt >() )
    {
      ic();             // set initial conditions
      advance();        // start time stepping
    }

    //! Set initial conditions
    void ic() {
      for (const auto& eq : g_diffeqs) eq.initialize( m_particles );
      m_proxy.init();   // signal to host that initialization is complete
    }

    //! Advance all equations one step in time
    void advance() {
      for (const auto& e : g_diffeqs) e.advance( m_particles, CkMyPe(), m_dt );
      ++m_it;
      m_t += m_dt;
      const auto term = g_inputdeck.get< tag::discr, tag::term >();
      if (m_t > term) m_t = term;
      if (std::fabs(m_t - term) > std::numeric_limits<tk::real>::epsilon() &&
          m_it < g_inputdeck.get< tag::discr, tag::nstep >())
        m_proxy.step();         // signal to host that time step is complete
      else
        m_proxy.finish();       // signal to host that time stepping is finished
    }

  private:
    Proxy m_proxy;              //!< Host proxy
    ParProps m_particles;       //!< Particle properties operated on
    uint64_t m_it;              //!< Iteration count
    tk::real m_t;               //!< Time
    tk::real m_dt;              //!< Time step size
};

} // quinoa::

#define CK_TEMPLATES_ONLY
#include <integrator.def.h>
#undef CK_TEMPLATES_ONLY

#endif // Integrator_h
