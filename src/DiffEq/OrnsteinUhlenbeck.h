// *****************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     System of Ornstein-Uhlenbeck SDEs
  \details   This file implements the time integration of a system of stochastic
    differential equations (SDEs), with linear drift and constant diffusion,
    whose invariant is the [joint normal
    distribution](http://en.wikipedia.org/wiki/Multivariate_normal_distribution).

    In a nutshell, the equation integrated governs a set of scalars,
    \f$Y_\alpha\f$, \f$\alpha\!=\!1,\dots,N\f$, as
    \f[
       \mathrm{d}Y_\alpha(t) = \theta_\alpha\left(\mu_\alpha - Y_\alpha\right)
       \mathrm{d}t + \sum_{\beta=1}^N \sigma_{\alpha\beta}\mathrm{d}W_\beta(t),
       \qquad \alpha=1,\dots,N
    \f]
    with parameter vectors \f$\theta_\alpha > 0\f$, \f$\mu_\alpha\f$, and
    symmetric positive semi-definite diffusion matrix
    \f$\sigma_{\alpha\beta}\f$. Here \f$\mathrm{d}W_\beta(t)\f$ is an isotropic
    vector-valued [Wiener process](http://en.wikipedia.org/wiki/Wiener_process)
    with independent increments. The invariant distribution is the joint normal
    distribution. This system of SDEs consists of N coupled equations, each
    being a single-variate [Ornstein-Uhlenbeck
    process](http://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process).

    From the Fokker-Planck equation, equivalent to the SDE above, the equations
    governing the means, \f$ \langle Y_\alpha \rangle\f$, are
    \f[
      \newcommand{\irmean}[1]{{\langle{#1}\rangle}}
      \frac{\partial\irmean{Y_\alpha}}{\partial t} =
        \theta_\alpha\left(\mu_\alpha - \irmean{Y_\alpha}\right)
    \f]
    while the equation governing the covariance matrix, \f$ \langle y_\alpha
    y_\beta \rangle \equiv \left\langle (Y_\alpha - \langle Y_\alpha \rangle)
    (Y_\beta - \langle Y_\beta\rangle) \right\rangle \f$, is
    \f[
      \newcommand{\irmean}[1]{{\langle{#1}\rangle}}
      \newcommand{\irv}[1]{\langle{#1^2}\rangle}
      \frac{\partial\irmean{y_\alpha y_\beta}}{\partial t} =
      -\left(\theta_\alpha+\theta_\beta\right)\irmean{y_\alpha y_\beta}
      +\sum_{\gamma=1}^N \sigma_{\alpha\gamma} \sigma_{\gamma\beta}.
    \f]
*/
// *****************************************************************************
#ifndef OrnsteinUhlenbeck_h
#define OrnsteinUhlenbeck_h

#include <vector>
#include <cmath>

#include "QuinoaConfig.h"

#ifdef HAS_MKL
  #include "NoWarning/mkl_lapacke.h"
#else
  #include "NoWarning/lapacke.h"
#endif

#include "InitPolicy.h"
#include "OrnsteinUhlenbeckCoeffPolicy.h"
#include "RNG.h"
#include "Particles.h"
#include "SystemComponents.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! \brief Ornstein-Uhlenbeck SDE used polymorphically with DiffEq
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Init - initialization policy, see DiffEq/InitPolicy.h
//!   - Coefficients - coefficients policy, see
//!       DiffEq/OrnsteinUhlenbeckCoeffPolicy.h
template< class Init, class Coefficients >
class OrnsteinUhlenbeck {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of Ornstein-Uhlenbeck SDEs to
    //!   construct. There can be multiple ornstein-uhlenbeck ... end blocks in
    //!   a control file. This index specifies which Ornstein-Uhlenbeck SDE
    //!   system to instantiate. The index corresponds to the order in which the
    //!   ornstein-uhlenbeck ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit OrnsteinUhlenbeck( ncomp_t c ) :
      m_c( c ),
      m_depvar( g_inputdeck.get< tag::param, tag::ou, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::ou >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::ou >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::ou, tag::rng >().at(c) ) ) ),
      m_sigma(),
      m_theta(),
      m_mu(),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::ou, tag::sigmasq >().at(c),
             g_inputdeck.get< tag::param, tag::ou, tag::theta >().at(c),
             g_inputdeck.get< tag::param, tag::ou, tag::mu >().at(c),
             m_sigma, m_theta, m_mu )
    {
      // Compute diffusion matrix using Cholesky-decomposition
      lapack_int n = static_cast< lapack_int >( m_ncomp );
      #ifndef NDEBUG
      lapack_int info =
      #endif
        LAPACKE_dpotrf( LAPACK_ROW_MAJOR, 'U', n, m_sigma.data(), n );
      Assert( info == 0, "Error in Cholesky-decomposition" );
    }

    //! Initalize SDE, prepare for time integration
    //! \param[in] stream Thread (or more precisely stream) ID 
    //! \param[in,out] particles Array of particle properties 
    //! \author J. Bakosi
    void initialize( int stream, tk::Particles& particles ) {
      //! Set initial conditions using initialization policy
      Init::template
        init< tag::ou >
            ( g_inputdeck, m_rng, stream, particles, m_c, m_ncomp, m_offset );
    }

    //! \brief Advance particles according to the system of Orsntein-Uhlenbeck
    //!   SDEs
    //! \param[in,out] particles Array of particle properties
    //! \param[in] stream Thread (or more precisely stream) ID
    //! \param[in] dt Time step size
    //! \author J. Bakosi
    void advance( tk::Particles& particles,
                  int stream,
                  tk::real dt,
                  tk::real,
                  const std::map< tk::ctr::Product, tk::real >& )
    {
      const auto npar = particles.nunk();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        std::vector< tk::real > dW( m_ncomp );
        m_rng.gaussian( stream, m_ncomp, dW.data() );

        // Advance all m_ncomp scalars
        for (ncomp_t i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          par += m_theta[i]*(m_mu[i] - par)*dt;
          for (ncomp_t j=0; j<m_ncomp; ++j) {
            tk::real d = m_sigma[ j*m_ncomp+i ] * sqrt(dt);     // use transpose
            par += d*dW[j];
          }
        }
      }
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const char m_depvar;                //!< Dependent variable
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    //! Coefficients
    std::vector< kw::sde_sigmasq::info::expect::type > m_sigma;
    std::vector< kw::sde_theta::info::expect::type > m_theta;
    std::vector< kw::sde_mu::info::expect::type > m_mu;

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // OrnsteinUhlenbeck_h
