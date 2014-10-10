//******************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck.h
  \author    J. Bakosi
  \date      Fri 10 Oct 2014 03:16:05 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Ornstein-Uhlenbeck SDE
  \details   Ornstein-Uhlenbeck SDE.
*/
//******************************************************************************
#ifndef OrnsteinUhlenbeck_h
#define OrnsteinUhlenbeck_h

#include <cmath>

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
      ErrChk( sigma.size() > c, "Wrong number of OU SDE parameters 'sigma'");
      ErrChk( theta.size() > c, "Wrong number of OU SDE parameters 'theta'");
      ErrChk( mu.size() > c, "Wrong number of OU SDE parameters 'mu'");
      // Use coefficients policy to initialize coefficients
      Coefficients( m_ncomp, sigma[c], theta[c], mu[c], m_sigma, m_theta, m_mu );
    }

    //! Set initial conditions
    void initialize( ParProps& particles ) const { Init( { particles } ); }

    //! Advance particles
    void advance( ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance all m_ncomp scalars
        for (unsigned int i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_sigma[i] * m_sigma[i] * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          par += m_theta[i]*(m_mu[i] - par)*dt + d*dW[i];
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
