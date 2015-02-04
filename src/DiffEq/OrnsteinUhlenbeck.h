//******************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:46:51 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
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
*/
//******************************************************************************
#ifndef OrnsteinUhlenbeck_h
#define OrnsteinUhlenbeck_h

#include <cmath>

#ifdef HAS_MKL
  #include <mkl_lapacke.h>
#else
  #include <lapacke.h>
#endif

#include <InitPolicy.h>
#include <OrnsteinUhlenbeckCoeffPolicy.h>
#include <RNG.h>

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

  public:
    //! \brief Constructor
    //! \param[in] c Index specifying which system of Ornstein-Uhlenbeck SDEs to
    //!   construct. There can be multiple ornstein-uhlenbeck ... end blocks in
    //!   a control file. This index specifies which Ornstein-Uhlenbeck SDE
    //!   system to instantiate. The index corresponds to the order in which the
    //!   ornstein-uhlenbeck ... end blocks are given the control file.
    //! \author J. Bakosi
    explicit OrnsteinUhlenbeck( unsigned int c ) :
      m_depvar( g_inputdeck.get< tag::param, tag::ou, tag::depvar >().at(c) ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::ou >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::ou >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::ou, tag::rng >().at(c) ) ) ),
      coeff( m_ncomp,
             g_inputdeck.get< tag::param, tag::ou, tag::sigma >().at(c),
             g_inputdeck.get< tag::param, tag::ou, tag::theta >().at(c),
             g_inputdeck.get< tag::param, tag::ou, tag::mu >().at(c),
             m_sigma, m_theta, m_mu )
    {
      // Compute diffusion matrix using Cholesky-decomposition
      lapack_int n = static_cast< lapack_int >( m_ncomp );
      lapack_int info =
        LAPACKE_dpotrf( LAPACK_ROW_MAJOR, 'U', n, m_sigma.data(), n );
      Assert( info == 0, "Error in Cholesky-decomposition" );
    }

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

    //! \brief Advance particles according to the system of Orsntein-Uhlenbeck
    //!   SDEs
    //! \author J. Bakosi
    void advance( tk::ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[ m_ncomp ];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance all m_ncomp scalars
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          par += m_theta[i]*(m_mu[i] - par)*dt;
          for (tk::ctr::ncomp_type j=0; j<m_ncomp; ++j) {
            tk::real d = m_sigma[ j*m_ncomp+i ] * sqrt(dt);     // use transpose
            par += d*dW[j];
          }
        }
      }
    }

  private:
    const char m_depvar;                //!< Dependent variable
    const tk::ctr::ncomp_type m_ncomp;  //!< Number of components
    const int m_offset;                 //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator

    //! Coefficients
    std::vector< kw::sde_sigma::info::expect::type > m_sigma;
    std::vector< kw::sde_theta::info::expect::type > m_theta;
    std::vector< kw::sde_mu::info::expect::type > m_mu;

    //! Coefficients policy
    Coefficients coeff;
};

} // walker::

#endif // OrnsteinUhlenbeck_h
