//******************************************************************************
/*!
  \file      src/SDE/GenDirichlet.h
  \author    J. Bakosi
  \date      Fri 10 Oct 2014 03:18:19 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Lochner's generalized Dirichlet SDE
  \details   Lochner's generalized Dirichlet SDE,
             see http://dx.doi.org/10.1063/1.4822416
*/
//******************************************************************************
#ifndef GenDirichlet_h
#define GenDirichlet_h

#include <InitPolicy.h>
#include <GenDirCoeffPolicy.h>
#include <RNG.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck;
extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! Lochner's generalized Dirichlet SDE used polymorphically with DiffEq
template< class Init, class Coefficients >
class GenDirichlet {

  public:
    //! Constructor: use coefficients policy to initialize coefficients
    explicit GenDirichlet( unsigned int c ) :
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::gendir >()[c] ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::gendir >(c) ),
      m_rng( g_rng.at( tk::ctr::raw(
        g_inputdeck.get< tag::param, tag::gendir, tk::tag::rng >()[c] ) ) )
    {
      const auto& b = g_inputdeck.get< tag::param, tag::gendir, tag::b >();
      const auto& S = g_inputdeck.get< tag::param, tag::gendir, tag::S >();
      const auto& k = g_inputdeck.get< tag::param, tag::gendir, tag::kappa >();
      const auto& C = g_inputdeck.get< tag::param, tag::gendir, tag::c >();
      ErrChk( b.size() > c, "Wrong number of OU SDE parameters 'b'");
      ErrChk( S.size() > c, "Wrong number of OU SDE parameters 'S'");
      ErrChk( k.size() > c, "Wrong number of OU SDE parameters 'kappa'");
      ErrChk( C.size() > c, "Wrong number of OU SDE parameters 'c'");
      // Use coefficients policy to initialize coefficients
      Coefficients( m_ncomp, b[c], S[c], k[c], C[c], m_b, m_S, m_k, m_c );
    }

    //! Set initial conditions
    void initialize( ParProps& particles ) const { Init( { particles } ); }

    //! Advance particles
    void advance( ParProps& particles, int stream, tk::real dt ) const {
      const auto npar = particles.npar();
      for (auto p=decltype(npar){0}; p<npar; ++p) {
        // Y_i = 1 - sum_{k=1}^{i} y_k
        tk::real Y[m_ncomp];
        Y[0] = 1.0 - particles( p, 0, m_offset );
        for (unsigned int i=1; i<m_ncomp; ++i)
          Y[i] = Y[i-1] - particles( p, i, m_offset );

        // U_i = prod_{j=1}^{K-i} 1/Y_{K-j}
        tk::real U[m_ncomp];
        U[m_ncomp-1] = 1.0;
        for (int i=m_ncomp-2; i>=0; --i) U[i] = U[i+1]/Y[i];

        // Generate Gaussian random numbers with zero mean and unit variance
        tk::real dW[m_ncomp];
        m_rng.gaussian( stream, m_ncomp, dW );

        // Advance first m_ncomp (K=N-1) scalars
        int k=0;
        for (unsigned int i=0; i<m_ncomp; ++i) {
          tk::real& par = particles( p, i, m_offset );
          tk::real d = m_k[i] * par * Y[m_ncomp-1] * U[i] * dt;
          d = (d > 0.0 ? sqrt(d) : 0.0);
          tk::real a=0.0;
          for (unsigned int j=i; j<m_ncomp-1; ++j) a += m_c[k++]/Y[j];
          par += U[i]/2.0*( m_b[i]*( m_S[i]*Y[m_ncomp-1] - (1.0-m_S[i])*par ) +
                            par*Y[m_ncomp-1]*a )*dt + d*dW[i];
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
    std::vector< tk::real > m_c;
};

} // quinoa::

#endif // GenDirichlet_h
