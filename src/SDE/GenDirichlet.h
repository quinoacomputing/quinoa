//******************************************************************************
/*!
  \file      src/SDE/GenDirichlet.h
  \author    J. Bakosi
  \date      Thu 06 Feb 2014 05:48:56 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Lochner's generalized Dirichlet SDE
  \details   Lochner's generalized Dirichlet SDE,
             see http://dx.doi.org/10.1063/1.4822416
*/
//******************************************************************************
#ifndef GenDirichlet_h
#define GenDirichlet_h

#include <SDE.h>
#include <GenDirCoeffPolicy.h>

namespace quinoa {

//! Lochner's generalized Dirichlet
template< class Init, class Coefficients >
class GenDirichlet : public SDE< Init, Coefficients > {

  public:
    //! SDE base shorthand
    using sde = SDE< Init, Coefficients >;

    //! Constructor
    explicit GenDirichlet( const Base& base,
                           const ParProps& particles,
                           int offset,
                           int ncomp ) :
      sde( base,
           base.control.get< tag::param, tag::gendir, tk::tag::rng >(),
           particles,
           offset,
           ncomp ),
      m_b( base.control.get< tag::component, tag::nscalar >() ),
      m_S( base.control.get< tag::component, tag::nscalar >() ),
      m_k( base.control.get< tag::component, tag::nscalar >() ),
      m_c( base.control.get< tag::component, tag::nscalar >() ),
      m_coeff( base.control.get< tag::component, tag::nscalar >(),
               base.control.get< tag::param, tag::gendir, tag::b >(),
               base.control.get< tag::param, tag::gendir, tag::S >(),
               base.control.get< tag::param, tag::gendir, tag::kappa >(),
               base.control.get< tag::param, tag::gendir, tag::c >(),
               m_b, m_S, m_k, m_c ) {}

    //! Pull base class data to scope
    using sde::m_particles;
    using sde::m_offset;
    using sde::m_ncomp;
    using sde::m_rng;

    //! Advance particles
    void advance(int p, int tid, tk::real dt) override {
      // Y_i = 1 - sum_{k=1}^{i} y_k
      tk::real Y[m_ncomp];
      Y[0] = 1.0 - m_particles( p, 0, m_offset );
      for (int i=1; i<m_ncomp; ++i) {
        Y[i] = Y[i-1] - m_particles( p, i, m_offset );
      }

      // U_i = prod_{j=1}^{K-i} 1/Y_{K-j}
      tk::real U[m_ncomp];
      U[m_ncomp-1] = 1.0;
      for (int i=m_ncomp-2; i>=0; --i) U[i] = U[i+1]/Y[i];

      // Generate Gaussian random numbers with zero mean and unit variance
      tk::real dW[m_ncomp];
      m_rng->gaussian( tid, m_ncomp, dW );

      // Advance first m_ncomp (K=N-1) scalars
      int k=0;
      for (int i=0; i<m_ncomp; ++i) {
        tk::real& par = m_particles( p, i, m_offset );
        tk::real d = m_k[i] * par * Y[m_ncomp-1] * U[i] * dt;
        d = (d > 0.0 ? sqrt(d) : 0.0);
        tk::real a=0.0;
        for (int j=i; j<m_ncomp-1; ++j) a += m_c[k++]/Y[j];
        par += U[i]/2.0*( m_b[i]*( m_S[i]*Y[m_ncomp-1] - (1.0-m_S[i])*par ) +
                          par*Y[m_ncomp-1]*a )*dt + d*dW[i];
      }
    }

  private:
    //! Don't permit copy constructor
    GenDirichlet(const GenDirichlet&) = delete;
    //! Don't permit copy assigment
    GenDirichlet& operator=(const GenDirichlet&) = delete;
    //! Don't permit move constructor
    GenDirichlet(GenDirichlet&&) = delete;
    //! Don't permit move assigment
    GenDirichlet& operator=(GenDirichlet&&) = delete;

    std::vector< tk::real > m_b;        //!< SDE coefficients
    std::vector< tk::real > m_S;
    std::vector< tk::real > m_k;
    std::vector< tk::real > m_c;

    Coefficients m_coeff;               //!< Coefficients policy
};

} // quinoa::

#endif // GenDirichlet_h
