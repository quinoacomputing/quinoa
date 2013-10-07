//******************************************************************************
/*!
  \file      src/Model/Mix/GenDirichlet/GenDirichlet.C
  \author    J. Bakosi
  \date      Mon Oct  7 10:34:44 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Mix/Mix.h>
#include <Mix/GenDirichlet/GenDirichlet.h>
#include <JPDF.h>

using namespace std;
using namespace quinoa;

void
GeneralizedDirichlet::advance(int p, int tid, const tk::real dt)
//******************************************************************************
//  Advance particles with the generalized Dirichlet model
//! \param[in]  p    Particle to advance
//! \param[in]  tid  Thread id
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  int i, j, k;
  tk::real d, a;
  tk::real* y;
  tk::real dW[m_nscalar];
  tk::real Y[m_nscalar];    //!< Y_i = 1 - sum_{k=1}^{i} y_k
  tk::real U[m_nscalar];    //!< U_i = prod_{j=1}^{K-i} 1/Y_{K-j}

  // Get access to particle scalars
  y = m_particles + p*m_nprop + m_offset;

  Y[0] = 1.0 - y[0];
  for (i=1; i<m_nscalar; ++i) Y[i] = Y[i-1] - y[i];

  U[m_nscalar-1] = 1.0;
  for (i=m_nscalar-2; i>=0; --i) U[i] = U[i+1]/Y[i];

  // Generate Gaussian random numbers with zero mean and unit variance
  rndstr()->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                     m_str[tid], m_nscalar, dW, 0.0, 1.0);

  // Advance first m_nscalar (K=N-1) scalars
  k=0;
  for (i=0; i<m_nscalar; ++i) {
    d = m_k[i]*y[i]*Y[m_nscalar-1]*U[i]*dt;
    if (d > 0.0) d = sqrt(d); else d = 0.0;
    a=0.0;
    for (j=i; j<m_nscalar-1; ++j) a += m_c[k++]/Y[j];
    y[i] += U[i]/2.0*(m_b[i]*(m_S[i]*Y[m_nscalar-1] - (1.0-m_S[i])*y[i]) +
                      y[i]*Y[m_nscalar-1]*a)*dt + d*dW[i];
  }
}

void
GeneralizedDirichlet::jpdf(JPDF& jpdf)
//******************************************************************************
//  Estimate joint scalar probability density function
//! \author  J. Bakosi
//******************************************************************************
{
IGNORE(jpdf);
//   for (uint64_t p=0; p<m_npar; ++p) {
//     tk::real* y = m_scalars + p*m_nscalar;
//     vector<tk::real> v(y, y+m_nscalar);
//     jpdf.insert(v);
//   }
}
