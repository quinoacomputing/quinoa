//******************************************************************************
/*!
  \file      src/DiffEq/Beta.h
  \author    J. Bakosi
  \date      Fri 13 Feb 2015 02:49:16 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     System of beta SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs) with linear drift and quadratic diagonal
    diffusion, whose invariant is the joint [beta
    distribution](http://en.wikipedia.org/wiki/Beta_distribution).

    In a nutshell, the equation integrated governs a set of scalars,
    \f$0\!\le\!Y_\alpha\f$, \f$\alpha\!=\!1,\dots,N\f$, as
    \f[
       \mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\left(S_\alpha - Y_\alpha\right)
       \mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha(1-Y_\alpha)}
       \mathrm{d}W_\alpha(t), \qquad \alpha=1,\dots,N
    \f]
    with parameter vectors \f$b_\alpha > 0\f$, \f$\kappa_\alpha > 0\f$, and \f$0
    < S_\alpha < 1\f$. Here
    \f$\mathrm{d}W_\alpha(t)\f$ is an isotropic vector-valued [Wiener
    process](http://en.wikipedia.org/wiki/Wiener_process) with independent
    increments. The invariant distribution is the joint beta distribution. This
    system of SDEs consists of N independent equations. For
    more on the beta SDE, see http://dx.doi.org/10.1080/14685248.2010.510843.
*/
//******************************************************************************
#ifndef Beta_h
#define Beta_h

#include <cmath>

#include <InitPolicy.h>
#include <BetaCoeffPolicy.h>
#include <RNG.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Beta SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see DiffEq/BetaCoeffPolicy.h
template< class Init, class Coefficients >
class Beta {

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of beta SDEs to construct.
    //!   There can be multiple beta ... end blocks in a control file. This
    //!   index specifies which beta SDE system to instantiate. The index
    //!   corresponds to the order in which the beta ... end blocks are given
    //!   the control file.
    //! \author J. Bakosi
    explicit Beta( unsigned int c ) :
      m_depvar( g_inputdeck.get< tag::param, tag::beta, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::beta >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::beta >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::beta, tag::rng >().at(c) ) ) ),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::beta, tag::b >().at(c),
             g_inputdeck.get< tag::param, tag::beta, tag::S >().at(c),
             g_inputdeck.get< tag::param, tag::beta, tag::kappa >().at(c),
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

    //! \brief Advance particles according to the system of beta SDEs
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
          tk::real d = m_k[i] * par * (1.0 - par) * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          par += 0.5*m_b[i]*(m_S[i] - par)*dt + d*dW[i];
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

#endif // Beta_h
