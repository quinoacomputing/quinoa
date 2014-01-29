//******************************************************************************
/*!
  \file      src/SDE/Dirichlet.h
  \author    J. Bakosi
  \date      Wed Jan 29 14:56:55 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet SDE
  \details   Dirichlet SDE, see http://dx.doi.org/10.1155/2013/842981
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <SDE.h>
#include <DirCoeffPolicy.h>

namespace quinoa {

//! Dirichlet
template< class Init, class Coefficients >
class Dirichlet : public SDE< Init > {

  public:
    //! SDE base shorthand
    using sde = SDE< Init >;

    //! Constructor
    explicit Dirichlet( const Base& base, const ParProps& particles ) :
      sde( base,
           base.control.get< tag::param, tag::dirichlet, tk::tag::rng >(),
           particles,
           base.control.scalarOffset(),
           base.control.get< tag::component, tag::nscalar >() ),
      m_b( base.control.get< tag::component, tag::nscalar >() ),
      m_S( base.control.get< tag::component, tag::nscalar >() ),
      m_k( base.control.get< tag::component, tag::nscalar >() ),
      m_coeff( m_coeffPolicy,
               base.control.get< tag::component, tag::nscalar >(),
               base.control.get< tag::param, tag::dirichlet, tag::b >(),
               base.control.get< tag::param, tag::dirichlet, tag::S >(),
               base.control.get< tag::param, tag::dirichlet, tag::kappa >(),
               m_b, m_S, m_k ) {}

    //! Return coefficients policy
    const std::string& coeffPolicy() const noexcept override {
      return m_coeffPolicy;
    }

    //! Pull base class data to scope
    using sde::m_particles;
    using sde::m_nprop;
    using sde::m_offset;
    using sde::m_ncomp;
    using sde::m_rng;

    //! Advance particles
    void advance(int p, int tid, tk::real dt) override {
      // Compute Nth scalar
      tk::real yn = 1.0 - m_particles(p, 0, m_offset);
      for (int i=1; i<m_ncomp; ++i) yn -= m_particles(p, i, m_offset);

      // Generate Gaussian random numbers with zero mean and unit variance
      tk::real dW[m_ncomp];
      m_rng->gaussian( tid, m_ncomp, dW );

      // Advance first m_ncomp (K=N-1) scalars
      for (int i=0; i<m_ncomp; ++i) {
        tk::real& par = m_particles( p, i, m_offset );
        tk::real d = m_k[i] * par * yn * dt;
        d = (d > 0.0 ? sqrt(d) : 0.0);
        par += 0.5*m_b[i]*( m_S[i]*yn - (1.0-m_S[i]) * par )*dt + d*dW[i];
      }
    }

  private:
    //! Don't permit copy constructor
    Dirichlet(const Dirichlet&) = delete;
    //! Don't permit copy assigment
    Dirichlet& operator=(const Dirichlet&) = delete;
    //! Don't permit move constructor
    Dirichlet(Dirichlet&&) = delete;
    //! Don't permit move assigment
    Dirichlet& operator=(Dirichlet&&) = delete;

    std::vector< tk::real > m_b;        //!< SDE coefficients
    std::vector< tk::real > m_S;
    std::vector< tk::real > m_k;

    std::string m_coeffPolicy;          //!< Coefficients policy name
    Coefficients m_coeff;               //!< Coefficients policy
};

} // quinoa::

#endif // Dirichlet_h
