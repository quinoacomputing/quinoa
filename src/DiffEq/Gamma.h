//******************************************************************************
/*!
  \file      src/DiffEq/Gamma.h
  \author    J. Bakosi
  \date      Fri 13 Feb 2015 02:49:47 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     System of gamma SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs), with linear drift and linear diagonal
    diffusion, whose invariant is the joint [gamma
    distribution](http://en.wikipedia.org/wiki/Gamma_distribution).

    In a nutshell, the equation integrated governs a set of scalars,
    \f$0\!\le\!Y_\alpha\f$, \f$\alpha\!=\!1,\dots,N\f$, as
    \f[
       \mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\big[S_\alpha -
       (1-S_\alpha)Y_\alpha\big]\mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha}
       \mathrm{d}W_\alpha(t), \qquad \alpha=1,\dots,N
    \f]
    with parameter vectors \f$b_\alpha > 0\f$, \f$\kappa_\alpha > 0\f$, and \f$0
    < S_\alpha < 1\f$. Here \f$\mathrm{d}W_\alpha(t)\f$ is an isotropic
    vector-valued [Wiener process](http://en.wikipedia.org/wiki/Wiener_process)
    with independent increments. The invariant distribution is the joint gamma
    distribution. This system of SDEs consists of N independent equations.
*/
//******************************************************************************
#ifndef Gamma_h
#define Gamma_h

#include <cmath>

#include <InitPolicy.h>
#include <GammaCoeffPolicy.h>
#include <RNG.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Gamma SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/GammaCoeffPolicy.h
template< class Init, class Coefficients >
class Gamma {

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of gamma SDEs to construct.
    //!   There can be multiple gamma ... end blocks in a control file. This
    //!   index specifies which gamma SDE system to instantiate. The index
    //!   corresponds to the order in which the gamma ... end blocks are given
    //!   the control file.
    //! \author J. Bakosi
    explicit Gamma( unsigned int c ) :
      m_depvar( g_inputdeck.get< tag::param, tag::gamma, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::gamma >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::gamma >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::gamma, tag::rng >().at(c) ) ) ),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::gamma, tag::b >().at(c),
             g_inputdeck.get< tag::param, tag::gamma, tag::S >().at(c),
             g_inputdeck.get< tag::param, tag::gamma, tag::kappa >().at(c),
             m_b, m_S, m_k ) {}

    //! Initalize SDE, prepare for time integration
    //! \param[inout] particles Array of particle properties 
    //! \param[in] stat Statistics object for accessing moments 
    //! \author J. Bakosi
    void initialize( tk::ParProps& particles, const tk::Statistics& stat ) {
      //! Set initial conditions using initialization policy
      Init( { particles } );
      //! Pre-lookup required statistical moments
      coeff.lookup( stat, m_depvar );
    }

    //! \brief Advance particles according to the system of gamma SDEs
    //! \author J. Bakosi
    void advance( tk::ParProps& particles, int stream, tk::real dt ) {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance all m_ncomp scalars
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_k[i] * par * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          par += 0.5*m_b[i]*(m_S[i] - (1.0 - m_S[i])*par)*dt + d*dW[i];
        }
      }
    }

  private:
    const char m_depvar;                //!< Dependent variable
    const tk::ctr::ncomp_type m_ncomp;  //!< Number of components
    const int m_offset;                 //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    //! Coefficients
    std::vector< kw::sde_b::info::expect::type > m_b;
    std::vector< kw::sde_S::info::expect::type > m_S;
    std::vector< kw::sde_kappa::info::expect::type > m_k;

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // Gamma_h
