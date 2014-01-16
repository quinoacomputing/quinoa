//******************************************************************************
/*!
  \file      src/SDE/Dirichlet.h
  \author    J. Bakosi
  \date      Wed 15 Jan 2014 09:44:34 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet SDE
  \details   Dirichlet SDE
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <SDE.h>
#include <DirCoeffPolicy.h>

namespace quinoa {

//! Dirichlet : Mix
template< class Init, class Layout, class Coefficients >
class Dirichlet : public SDE< Init, Layout > {

  public:
    //! SDE base shorthand
    using sde = SDE< Init, Layout >;

    //! Constructor
    explicit Dirichlet( const Base& base, tk::real* const particles ) :
      sde( base,
           particles,
           base.control.scalarOffset(),
           base.control.get< ctr::component, ctr::nscalar >() ),
      m_coeff( base.control.get< ctr::component, ctr::nscalar >(),
               base.control.get< ctr::param, ctr::dirichlet, ctr::b >(),
               base.control.get< ctr::param, ctr::dirichlet, ctr::S >(),
               base.control.get< ctr::param, ctr::dirichlet, ctr::kappa >(),
               m_b, m_S, m_k ) {}

    //! Pull base class data to scope
    using sde::m_particles;
    using sde::m_nprop;
    using sde::m_offset;
    using sde::m_ncomp;

    //! Advance particles
    void advance(int p, int tid, tk::real dt) override {
      // Get access to particle scalars
      tk::real* y = m_particles + p*m_nprop + m_offset;

      // Compute Nth scalar
      tk::real yn = 1.0 - y[0];
      for (int i=1; i<m_ncomp; ++i) yn -= y[i];

      // Generate Gaussian random numbers with zero mean and unit variance
      tk::real dW[m_ncomp];
      //rndstr()->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
      //                    m_str[tid], m_ncomp, dW, 0.0, 1.0);

      // Advance first m_ncomp (K=N-1) scalars
      for (int i=0; i<m_ncomp; ++i) {
        tk::real d = m_k[i]*y[i]*yn*dt;
        if (d > 0.0) d = sqrt(d); else d = 0.0;
        y[i] += 0.5*m_b[i]*(m_S[i]*yn - (1.0-m_S[i])*y[i])*dt + d*dW[i];
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

    Coefficients m_coeff;               //!< Coefficients update policy

    std::vector< tk::real > m_b;        //!< SDE coefficients
    std::vector< tk::real > m_S;
    std::vector< tk::real > m_k;
};

} // quinoa::

#endif // Dirichlet_h
