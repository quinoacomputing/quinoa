//******************************************************************************
/*!
  \file      src/DiffEq/Gamma.h
  \author    J. Bakosi
  \date      Fri 05 Dec 2014 03:11:39 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Gamma SDE
  \details   Gamma SDE.
*/
//******************************************************************************
#ifndef Gamma_h
#define Gamma_h

#include <cmath>

#include <InitPolicy.h>
#include <GammaCoeffPolicy.h>
#include <RNG.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! Gamma SDE used polymorphically with DiffEq
template< class Init, class Coefficients >
class Gamma {

  public:
    //! Constructor
    explicit Gamma( unsigned int c ) :
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::gamma >()[c] ),
      m_offset(g_inputdeck.get< tag::component >().offset< tag::gamma >(c)),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::gamma, tk::tag::rng >()[c] ) ) )
    {
      const auto& b = g_inputdeck.get< tag::param, tag::gamma, tag::b >();
      const auto& S = g_inputdeck.get< tag::param, tag::gamma, tag::S >();
      const auto& k = g_inputdeck.get< tag::param, tag::gamma, tag::kappa >();
      ErrChk( b.size() > c, "Indexing out of gamma SDE parameters 'b'");
      ErrChk( S.size() > c, "Indexing out of gamma SDE parameters 'S'");
      ErrChk( k.size() > c, "Indexing out of gamma SDE parameters 'kappa'");
      // Use coefficients policy to initialize coefficients
      Coefficients( m_ncomp, b[c], S[c], k[c], m_b, m_S, m_k );
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
          tk::real d = m_k[i] * par * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          par += 0.5*m_b[i]*(m_S[i] - (1.0 - m_S[i])*par)*dt + d*dW[i];
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

} // quinoa::

#endif // Gamma_h
