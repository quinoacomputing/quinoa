//******************************************************************************
/*!
  \file      src/Model/Mass/Beta/Beta.C
  \author    J. Bakosi
  \date      Mon Oct 28 07:27:28 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Beta mass model
  \details   Beta mass model
*/
//******************************************************************************

#include <Mass/Beta/Beta.h>

// void
// Beta::init()
// //******************************************************************************
// //  Initialize densities
// //! \author  J. Bakosi
// //******************************************************************************
// {
//   //const tk::real a = 0.075;    // At = 0.5
//   //const tk::real a = 0.06;     // At = 0.25
//   //const tk::real a = 0.039;    // At = 0.05
// 
//   //initBeta(a, a, 0.0, 1.0);
//   //U[pP+9] = 2.0*At*x[0] + 1.0-At;
// }
// 
// void
// Beta::advance(int p, int tid, tk::real dt)
// //******************************************************************************
// //  Advance particles with the Beta model
// //! \param[in]  p    Particle to advance
// //! \param[in]  tid  Thread id
// //! \param[in]  dt   Time step size
// //! \author  J. Bakosi
// //******************************************************************************
// {
// IGNORE(p);
// IGNORE(tid);
// IGNORE(dt);
// }
// 
// void
// Beta::jpdf(tk::JPDF& jpdf)
// //******************************************************************************
// //  Estimate joint scalar probability density function
// //! \author  J. Bakosi
// //******************************************************************************
// {
// IGNORE(jpdf);
// //   for (int p=0; p<m_npar; ++p) {
// //     tk::real* y = m_scalars + p*m_nscalar;
// //     vector<tk::real> v(y, y+m_nscalar);
// //     jpdf.insert(v);
// //   }
// }
