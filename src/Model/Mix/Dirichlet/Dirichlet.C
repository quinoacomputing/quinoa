//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.C
  \author    J. Bakosi
  \date      Sun 12 May 2013 03:33:33 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Control.h>
#include <Mix.h>
#include <Dirichlet.h>
#include <JPDF.h>

using namespace std;
using namespace Quinoa;

Dirichlet::Dirichlet(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control,
                     real* const scalars) :
  Mix<Dirichlet>(memory,
                 paradigm,
                 control,
                 control->get<control::NSCALAR>(),
                 scalars),
  m_b(control->get<control::B>()),
  m_S(control->get<control::S>()),
  m_k(control->get<control::KAPPA>())
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
         "Wrong number of Dirichlet model parameters 'b'");
  ErrChk(m_S.size() == static_cast<unsigned int>(m_nscalar), FATAL,
         "Wrong number of Dirichlet model parameters 'S'");
  ErrChk(m_k.size() == static_cast<unsigned int>(m_nscalar), FATAL,
         "Wrong number of Dirichlet model parameters 'k'");
}

void
Dirichlet::echo() const
//******************************************************************************
//  Echo information on Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Dirichlet" << endl;
}

void
Dirichlet::init()
//******************************************************************************
//  Initialize scalars
//! \author  J. Bakosi
//******************************************************************************
{
  initUniform();
}

void
Dirichlet::advance(const real& dt)
//******************************************************************************
//  Advance particles with the Dirichlet model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  int tid, p, i;
  real yn, d;
  real* y;
  real dW[m_nscalar];

  #ifdef _OPENMP
  #pragma omp parallel private(tid, p, i, y, yn, dW, d)
  #endif
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

      // Compute Nth scalar
      yn = 1.0 - y[0];
      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (i=1; i<m_nscalar; ++i) yn -= y[i];

      // Generate Gaussian random numbers with zero mean and unit variance
      rndstr()->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                          m_str[tid], m_nscalar, dW, 0.0, 1.0);

      // Advance first m_nscalar (K=N-1) scalars
      for (i=0; i<m_nscalar; ++i) {
        d = m_k[i]*y[i]*yn*dt;
        if (d > 0.0) d = sqrt(d); else d = 0.0;
        y[i] += m_b[i]/2.0*(m_S[i]*yn - (1.0-m_S[i])*y[i])*dt + d*dW[i];
      }
    } // m_npar
  } // omp parallel
}

void
Dirichlet::jpdf(JPDF& jpdf)
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
