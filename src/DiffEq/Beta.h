//******************************************************************************
/*!
  \file      src/DiffEq/Beta.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 06:29:42 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Beta SDE
  \details   Beta SDE, see http://dx.doi.org/10.1080/14685248.2010.510843
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

//! Beta SDE used polymorphically with DiffEq
template< class Init, class Coefficients >
class Beta {

  public:
    //! Constructor
    explicit Beta( unsigned int c ) :
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::beta >()[c] ),
      m_offset(g_inputdeck.get< tag::component >().offset< tag::beta >(c)),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::beta, tag::rng >()[c] ) ) )
    {
      const auto& b = g_inputdeck.get< tag::param, tag::beta, tag::b >();
      const auto& S = g_inputdeck.get< tag::param, tag::beta, tag::S >();
      const auto& k = g_inputdeck.get< tag::param, tag::beta, tag::kappa >();
      ErrChk( b.size() > c, "Indexing out of beta SDE parameters 'b'");
      ErrChk( S.size() > c, "Indexing out of beta OU SDE parameters 'S'");
      ErrChk( k.size() > c, "Indexing out of beta OU SDE parameters 'kappa'");
      // Use coefficients policy to initialize coefficients
      Coefficients( m_ncomp, b[c], S[c], k[c], m_b, m_S, m_k );
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
        for (unsigned int i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_k[i] * par * (1.0 - par) * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          par += 0.5*m_b[i]*(m_S[i] - par)*dt + d*dW[i];
        }
      }
    }

  private:
    const unsigned int m_ncomp;         //!< Number of components
    const int m_offset;                 //!< Offset SDE operates from
    const tk::RNG& m_rng;               //!< Random number generator
    std::vector< tk::real > m_b;        //!< Coefficients
    std::vector< tk::real > m_S;
    std::vector< tk::real > m_k;
};

} // walker::

#endif // Beta_h
