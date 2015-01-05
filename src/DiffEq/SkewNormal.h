//******************************************************************************
/*!
  \file      src/DiffEq/SkewNormal.h
  \author    J. Bakosi
  \date      Tue 13 Jan 2015 10:53:33 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Skew-normal SDE
  \details   Skew-normal SDE.
*/
//******************************************************************************
#ifndef SkewNormal_h
#define SkewNormal_h

#include <cmath>

#include <InitPolicy.h>
#include <SkewNormalCoeffPolicy.h>
#include <RNG.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! Skew-normal SDE used polymorphically with DiffEq
template< class Init, class Coefficients >
class SkewNormal {

  public:
    //! Constructor
    explicit SkewNormal( unsigned int c ) :
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::skewnormal >()[c] ),
      m_offset(g_inputdeck.get< tag::component >().offset< tag::skewnormal >(c)),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::skewnormal, tag::rng >()[c] ) ) )
    {
      const auto& timescale =
        g_inputdeck.get< tag::param, tag::skewnormal, tag::timescale >();
      const auto& sigma =
        g_inputdeck.get< tag::param, tag::skewnormal, tag::sigma >();
      const auto& lambda =
        g_inputdeck.get< tag::param, tag::skewnormal, tag::lambda >();
      ErrChk( timescale.size() > c,
              "Indexing out of Skew-normal SDE parameters 'timescale'");
      ErrChk( sigma.size() > c,
              "Indexing out of Skew-normal SDE parameters 'sigma'");
      ErrChk( lambda.size() > c,
              "Indexing out of Skew-normal SDE parameters 'lambda'");
      // Use coefficients policy to initialize coefficients
      Coefficients( m_ncomp, timescale[c], sigma[c], lambda[c],
                    m_timescale, m_sigma, m_lambda );
    }

    //! Set initial conditions
    void initialize( tk::ParProps& particles ) const { Init( { particles } ); }

    //! Advance particles
    void advance( tk::ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance all m_ncomp scalars
        for (tk::ctr::ncomp_type i=0; i<m_ncomp; ++i) {
          tk::real& x = particles( p, i, m_offset );
          tk::real d = 2.0 * m_sigma[i] * m_sigma[i] / m_timescale[i] * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          x += - ( x - m_lambda[i] * m_sigma[i] * m_sigma[i]
                       * std::sqrt( 2.0 / pi() )
                       * std::exp( - m_lambda[i] * m_lambda[i] * x * x / 2.0 )
                       / ( 1.0 + std::erf( m_lambda[i] * x / std::sqrt(2.0) ) )
                   ) / m_timescale[i] * dt
                 + d*dW[i];
        }
      }
    }

  private:
    const tk::ctr::ncomp_type m_ncomp;    //!< Number of components
    const int m_offset;                   //!< Offset SDE operates from
    const tk::RNG& m_rng;                 //!< Random number generator
    std::vector< tk::real > m_timescale;  //!< Coefficients
    std::vector< tk::real > m_sigma;
    std::vector< tk::real > m_lambda;

    constexpr double pi() const { return std::atan(1.0) * 4.0; }
};

} // walker::

#endif // SkewNormal_h
