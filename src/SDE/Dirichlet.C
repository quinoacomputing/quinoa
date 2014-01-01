//******************************************************************************
/*!
  \file      src/SDE/Dirichlet.C
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:17:15 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet SDE
  \details   Dirichlet SDE
*/
//******************************************************************************

#include <Dirichlet.h>

using quinoa::Dirichlet;

// void
// Dirichlet::advance(int p, int tid, tk::real dt)
// //******************************************************************************
// //  Advance particles with the Dirichlet model
// //! \param[in]  p    Particle to advance
// //! \param[in]  tid  Thread id
// //! \param[in]  dt   Time step size
// //! \author  J. Bakosi
// //******************************************************************************
// {
//   int i;
//   tk::real yn, d;
//   tk::real* y;
//   tk::real dW[m_nscalar];
// 
//   // Get access to particle scalars
//   y = m_particles + p*m_nprop + m_offset;
// 
//   // Compute Nth scalar
//   yn = 1.0 - y[0];
//   for (i=1; i<m_nscalar; ++i) yn -= y[i];
// 
//   // Generate Gaussian random numbers with zero mean and unit variance
//   //rndstr()->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
//   //                    m_str[tid], m_nscalar, dW, 0.0, 1.0);
// 
//   // Advance first m_nscalar (K=N-1) scalars
//   for (i=0; i<m_nscalar; ++i) {
//     d = m_k[i]*y[i]*yn*dt;
//     if (d > 0.0) d = sqrt(d); else d = 0.0;
//     y[i] += 0.5*m_b[i]*(m_S[i]*yn - (1.0-m_S[i])*y[i])*dt + d*dW[i];
//   }
// }
// 
// void
// Dirichlet::jpdf(JPDF& jpdf)
// //******************************************************************************
// //  Estimate joint scalar probability density function
// //! \author  J. Bakosi
// //******************************************************************************
// {
// IGNORE(jpdf);
// //   for (uint64_t p=0; p<m_npar; ++p) {
// //     tk::real* y = m_scalars + p*m_nscalar;
// //     vector<tk::real> v(y, y+m_nscalar);
// //     jpdf.insert(v);
// //   }
// }
