//******************************************************************************
/*!
  \file      src/SDE/Dirichlet.h
  \author    J. Bakosi
  \date      Tue Jan 14 09:13:45 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet SDE
  \details   Dirichlet SDE
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <SDE.h>

namespace quinoa {

//! Dirichlet : Mix
template< class Init, class Coefficients >
class Dirichlet : public SDE< Init > {

  public:
    //! Constructor
    explicit Dirichlet(const Base& base, tk::real* const particles) :
      SDE< Init >
         ( base,
           particles,
           base.control.scalarOffset(),
           base.control.get<ctr::component, ctr::nscalar>() ) {}
//       m_b(base.control.get<ctr::param, ctr::dirichlet, ctr::b>()),
//       m_S(base.control.get<ctr::param, ctr::dirichlet, ctr::S>()),
//       m_k(base.control.get<ctr::param, ctr::dirichlet, ctr::kappa>()) {
//       ErrChk(m_b.size() == static_cast<unsigned int>(m_nscalar),
//              tk::ExceptType::FATAL,
//              "Wrong number of Dirichlet model parameters 'b'");
//       ErrChk(m_S.size() == static_cast<unsigned int>(m_nscalar),
//              tk::ExceptType::FATAL,
//              "Wrong number of Dirichlet model parameters 'S'");
//       ErrChk(m_k.size() == static_cast<unsigned int>(m_nscalar),
//              tk::ExceptType::FATAL,
//              "Wrong number of Dirichlet model parameters 'k'");
//     }

    //! Advance particles
    void advance(int p, int tid, tk::real dt) override {
//       int i;
//       tk::real yn, d;
//       tk::real* y;
//       tk::real dW[m_nscalar];
// 
//       // Get access to particle scalars
//       y = m_particles + p*m_nprop + m_offset;
// 
//       // Compute Nth scalar
//       yn = 1.0 - y[0];
//       for (i=1; i<m_nscalar; ++i) yn -= y[i];
// 
//       // Generate Gaussian random numbers with zero mean and unit variance
//       //rndstr()->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
//       //                    m_str[tid], m_nscalar, dW, 0.0, 1.0);
// 
//       // Advance first m_nscalar (K=N-1) scalars
//       for (i=0; i<m_nscalar; ++i) {
//         d = m_k[i]*y[i]*yn*dt;
//         if (d > 0.0) d = sqrt(d); else d = 0.0;
//         y[i] += 0.5*m_b[i]*(m_S[i]*yn - (1.0-m_S[i])*y[i])*dt + d*dW[i];
//       }
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

//     const std::vector<tk::real> m_b;         //!< SDE coefficients
//     const std::vector<tk::real> m_S;
//     const std::vector<tk::real> m_k;
};

} // quinoa::

#endif // Dirichlet_h
