//******************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck.h
  \author    J. Bakosi
  \date      Fri 05 Dec 2014 02:16:54 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Ornstein-Uhlenbeck SDE
  \details   Ornstein-Uhlenbeck SDE.
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
#include <OUCoeffPolicy.h>
#include <RNG.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! Ornstein-Uhlenbeck SDE used polymorphically with DiffEq
template< class Init, class Coefficients >
class OrnsteinUhlenbeck {

  public:
    //! Constructor
    explicit OrnsteinUhlenbeck( unsigned int c ) :
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::ou >()[c] ),
      m_offset(g_inputdeck.get< tag::component >().offset< tag::ou >(c)),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::ou, tk::tag::rng >()[c] ) ) )
    {
      const auto& sigma = g_inputdeck.get< tag::param, tag::ou, tag::sigma >();
      const auto& theta = g_inputdeck.get< tag::param, tag::ou, tag::theta >();
      const auto& mu = g_inputdeck.get< tag::param, tag::ou, tag::mu >();
      ErrChk( sigma.size() > c, "Indexing out of OU SDE parameters 'sigma'");
      ErrChk( theta.size() > c, "Indexing out of OU SDE parameters 'theta'");
      ErrChk( mu.size() > c, "Indexing out of OU SDE parameters 'mu'");
      // Use coefficients policy to initialize coefficients
      Coefficients( m_ncomp, sigma[c], theta[c], mu[c], m_sigma, m_theta, m_mu );

      // Compute diffusion matrix using Cholesky-decomposition
      lapack_int n = static_cast< lapack_int >( m_ncomp );
      lapack_int info =
        LAPACKE_dpotrf( LAPACK_ROW_MAJOR, 'U', n, m_sigma.data(), n );
    }

    //! Set initial conditions
    void initialize( ParProps& particles ) const { Init( { particles } ); }

    //! Advance particles
    void advance( ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        tk::real dW[ m_ncomp ];
        m_rng.gaussian( stream, m_ncomp, dW );
        // Advance all m_ncomp scalars
        for (unsigned int i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          par += m_theta[i]*(m_mu[i] - par)*dt;
          for (unsigned int j=0; j<m_ncomp; ++j) {
            tk::real d = m_sigma[ j*m_ncomp+i ] * sqrt(dt);     // use transpose
            par += d*dW[j];
          }
        }
      }
    }

  private:
    const unsigned int m_ncomp;         //!< Number of components
    const int m_offset;                 //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator
    std::vector< tk::real > m_sigma;    //!< Coefficients
    std::vector< tk::real > m_theta;
    std::vector< tk::real > m_mu;
};

} // quinoa::

#endif // OrnsteinUhlenbeck_h
