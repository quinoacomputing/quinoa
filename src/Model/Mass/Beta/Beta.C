//******************************************************************************
/*!
  \file      src/Model/Mass/Beta/Beta.C
  \author    J. Bakosi
  \date      Mon 13 May 2013 09:23:25 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Beta mass model
  \details   Beta mass model
*/
//******************************************************************************

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Control.h>
#include <Mass.h>
#include <Beta.h>
#include <JPDF.h>

using namespace std;
using namespace Quinoa;

Beta::Beta(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control,
                     real* const scalars) :
  Mass<Beta>(memory,
             paradigm,
             control,
             control->get<control::NSCALAR>(),
             scalars),
  m_At(control->get<control::AT>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  scalars  Pointer to particle scalars
//! \author  J. Bakosi
//******************************************************************************
{
}

void
Beta::init()
//******************************************************************************
//  Initialize scalars
//! \author  J. Bakosi
//******************************************************************************
{
  const real a = 0.075;    // At = 0.5
  //const real a = 0.06;     // At = 0.25
  //const real a = 0.039;    // At = 0.05

  initBeta(a, a, 0.0, 1.0);
  //U[pP+9] = 2.0*At*x[0] + 1.0-At;
  cout << m_At << endl;
}

void
Beta::advance(const real& dt)
//******************************************************************************
//  Advance particles with the Beta model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
//   int tid, p, i;
//   real yn, d;
//   real* y;
//   real dW[m_nscalar];
// 
//   #ifdef _OPENMP
//   #pragma omp parallel private(tid, p, i, y, yn, dW, d)
//   #endif
//   {
//     #ifdef _OPENMP
//     tid = omp_get_thread_num();
//     #else
//     tid = 0;
//     #endif
// 
//     #ifdef _OPENMP
//     #pragma omp for
//     #endif
//     for (p=0; p<m_npar; ++p) {
//       // Get access to particle scalars
//       y = m_scalars + p*m_nscalar;
// 
//       // Compute Nth scalar
//       yn = 1.0 - y[0];
//       #ifdef __INTEL_COMPILER
//       #pragma vector always
//       #endif
//       for (i=1; i<m_nscalar; ++i) yn -= y[i];
// 
//       // Generate Gaussian random numbers with zero mean and unit variance
//       rndstr()->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
//                           m_str[tid], m_nscalar, dW, 0.0, 1.0);
// 
//       // Advance first m_nscalar (K=N-1) scalars
//       for (i=0; i<m_nscalar; ++i) {
//         d = m_k[i]*y[i]*yn*dt;
//         if (d > 0.0) d = sqrt(d); else d = 0.0;
//         y[i] += m_b[i]/2.0*(m_S[i]*yn - (1.0-m_S[i])*y[i])*dt + d*dW[i];
//       }
//     } // m_npar
//   } // omp parallel
}

void
Beta::jpdf(JPDF& jpdf)
//******************************************************************************
//  Estimate joint scalar probability density function
//! \author  J. Bakosi
//******************************************************************************
{
//   for (int p=0; p<m_npar; ++p) {
//     real* y = m_scalars + p*m_nscalar;
//     vector<real> v(y, y+m_nscalar);
//     jpdf.insert(v);
//   }
}
