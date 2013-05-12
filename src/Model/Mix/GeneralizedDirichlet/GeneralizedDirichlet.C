//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.C
  \author    J. Bakosi
  \date      Sun 12 May 2013 03:41:11 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Control.h>
#include <Mix.h>
#include <GeneralizedDirichlet.h>
#include <JPDF.h>

using namespace std;
using namespace Quinoa;

GeneralizedDirichlet::GeneralizedDirichlet(Memory* const memory,
                                           Paradigm* const paradigm,
                                           Control* const control,
                                           real* const scalars) :
  Mix<GeneralizedDirichlet>(memory,
                            paradigm,
                            control,
                            control->get<control::NSCALAR>(),
                            scalars),
  m_b(control->get<control::B>()),
  m_S(control->get<control::S>()),
  m_k(control->get<control::KAPPA>()),
  m_c(control->get<control::C>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  scalars  Pointer to particle scalars
//! \author  J. Bakosi
//******************************************************************************
{
  ErrChk(m_b.size() == static_cast<unsigned int>(m_nscalar), FATAL,
            "Wrong number of generalized Dirichlet model parameters 'b'");
  ErrChk(m_S.size() == static_cast<unsigned int>(m_nscalar), FATAL, 
            "Wrong number of generalized Dirichlet model parameters 'S'");
  ErrChk(m_k.size() == static_cast<unsigned int>(m_nscalar), FATAL,
            "Wrong number of generalized Dirichlet model parameters 'k'");
//  ErrChk(m_c.size() == static_cast<unsigned int>(m_nscalar*(m_nscalar-1)/2),
//          FATAL, "Wrong number of generalized Dirichlet model parameters 'c'");
}

void
GeneralizedDirichlet::init()
//******************************************************************************
//  Initialize scalars
//! \author  J. Bakosi
//******************************************************************************
{
  initUniform();
}

void
GeneralizedDirichlet::advance(const real& dt)
//******************************************************************************
//  Advance particles with the generalized Dirichlet model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  int tid, p, i, j, k;
  real d, a;
  real* y;
  real dW[m_nscalar];
  real Y[m_nscalar];    //!< Y_i = 1 - sum_{k=1}^{i} y_k
  real U[m_nscalar];    //!< U_i = prod_{j=1}^{m_nscalar-i} 1/sY_{m_nscalar-j}

  #ifdef _OPENMP
  #pragma omp parallel private(tid, p, i, j, k, d, a, y, Y, U, dW)
  #endif // _OPENMP
  {
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<m_npar; ++p) {
      // Get access to particle scalars
      y = m_scalars + p*m_nscalar;

      Y[0] = 1.0 - y[0];
      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (i=1; i<m_nscalar; ++i) Y[i] = Y[i-1] - y[i];

      U[m_nscalar-1] = 1.0;
      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (i=m_nscalar-2; i>=0; --i) U[i] = U[i+1]/Y[i];

      // Generate Gaussian random numbers with zero mean and unit variance
      rndstr()->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                         m_str[tid], m_nscalar, dW, 0.0, 1.0);

      // Advance first m_nscalar (K=N-1) scalars
      k=0;
      a=0.0;
      for (i=0; i<m_nscalar; ++i) {
        d = m_k[i]*y[i]*Y[m_nscalar-1]*U[i]*dt;
        if (d > 0.0) d = sqrt(d); else d = 0.0;
        for (j=i; j<m_nscalar-1; ++j) a += m_c[k++]/Y[j];
        y[i] += U[i]/2.0*m_b[i]*
                ((m_S[i]*Y[m_nscalar-1] - (1.0-m_S[i])*y[i]) +
                  y[i]*Y[m_nscalar-1]*a)*dt + d*dW[i];
      }
    } // m_npar
  } // omp parallel
}

void
GeneralizedDirichlet::jpdf(JPDF& jpdf)
//******************************************************************************
//  Estimate joint scalar probability density function
//! \author  J. Bakosi
//******************************************************************************
{
  for (int p=0; p<m_npar; ++p) {
    real* y = m_scalars + p*m_nscalar;
    vector<real> v(y, y+m_nscalar);
    jpdf.insert(v);
  }
}
