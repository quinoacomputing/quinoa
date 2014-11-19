//******************************************************************************
/*!
  \file      src/DiffEq/Dirichlet.h
  \author    J. Bakosi
  \date      Mon 17 Nov 2014 07:46:18 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Dirichlet SDE
  \details   Dirichlet SDE, see http://dx.doi.org/10.1155/2013/842981,
             \f[\mathrm{d}Y_\alpha(t) = \frac{b_\alpha}{2}\big[S_\alpha Y_N -
             (1-S_\alpha)Y_\alpha\big]\mathrm{d}t + \sqrt{\kappa_\alpha Y_\alpha
             Y_N}\mathrm{d}W_\alpha(t), \qquad \alpha=1,\dots,N-1 \f]
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <cmath>

#include <InitPolicy.h>
#include <DirCoeffPolicy.h>
#include <RNG.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! Dirichlet SDE used polymorphically with DiffEq
template< class Init, class Coefficients >
class Dirichlet {

  public:
    //! Constructor
    explicit Dirichlet( unsigned int c ) :
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::dirichlet >()[c] ),
      m_offset(g_inputdeck.get< tag::component >().offset< tag::dirichlet >(c)),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::dirichlet, tk::tag::rng >()[c] ) ) )
    {
      const auto& b = g_inputdeck.get< tag::param, tag::dirichlet, tag::b >();
      const auto& S = g_inputdeck.get< tag::param, tag::dirichlet, tag::S >();
      const auto& k = g_inputdeck.get< tag::param, tag::dirichlet, tag::kappa >();
      ErrChk( b.size() > c, "Wrong number of Dirichlet SDE parameters 'b'");
      ErrChk( S.size() > c, "Wrong number of Dirichlet SDE parameters 'S'");
      ErrChk( k.size() > c, "Wrong number of Dirichlet SDE parameters 'kappa'");
      // Use coefficients policy to initialize coefficients
      Coefficients( m_ncomp, b[c], S[c], k[c], m_b, m_S, m_k );
    }

    //! Set initial conditions
    void initialize( ParProps& particles ) const { Init( { particles } ); }

    //! Advance particles
    void advance( ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Compute Nth scalar
        tk::real yn = 1.0 - particles(p, 0, m_offset);
        for (unsigned int i=1; i<m_ncomp; ++i) yn -= particles( p, i, m_offset );

        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance first m_ncomp (K=N-1) scalars
        for (unsigned int i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_k[i] * par * yn * dt;
          d = (d > 0.0 ? std::sqrt(d) : 0.0);
          par += 0.5*m_b[i]*( m_S[i]*yn - (1.0-m_S[i]) * par )*dt + d*dW[i];
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

#endif // Dirichlet_h
