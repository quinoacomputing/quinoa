//******************************************************************************
/*!
  \file      src/Model/Hydro/SLM/SLM.C
  \author    J. Bakosi
  \date      Sun 01 Sep 2013 03:08:23 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************

#include <iostream>
#include <cstring>
#include <cmath>

#include <Macro.h>
#include <Control.h>
#include <Hydro.h>
#include <SLM.h>

using namespace std;
using namespace quinoa;

void
SimplifiedLangevin::init()
//******************************************************************************
//  Initialize particle properties
//! \author  J. Bakosi
//******************************************************************************
{
//  initGaussian(m_velocities, m_nvelocity, rndstr(), m_str[0], 0.0, 1.0);
//   real r[m_nvelocity];
// 
//   // Generate initial values for all scalars for all particles
//   for (uint64_t p=0; p<m_npar; ++p) {
//     m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
//                        m_str[0], m_nvelocity, r, 0.0, 1.0);
//     memcpy(m_velocities + p*m_nvelocity, r, m_nvelocity*sizeof(real));
//   }
}

void
SimplifiedLangevin::advance(const real& dt)
//******************************************************************************
//  Advance particles with the simplified Langevin model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
IGNORE(dt);
//   uint64_t p;
//   int tid, i;
//   real d, tke, eps, S=1.0;
//   real* X;
//   real* U;
//   real dW[3], rs[6];
// 
//   #ifdef _OPENMP
//   //#pragma omp parallel private(tid, p, i, y, yn, dW, d)
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
//       // Get access to particle position, velocity
//       X = m_velocities + p*m_nvelocity;
//       U = X + 3;
// 
//       // Generate Gaussian random numbers with zero mean and unit variance
//       m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
//                          m_str[tid], 3, dW, 0.0, 1.0);
// 
//       tke = 0.5*(rs[0] + rs[1] + rs[2]);
//       eps = S*tke/2.36;
// 
//       // Advance velocity
//       d = m_C0*eps*dt;
//       if (d > 0.0) d = sqrt(d); else d = 0.0;
//       //U[0] += ()*dt + d*dW[i];
//     } // m_npar
//   } // omp parallel
}
