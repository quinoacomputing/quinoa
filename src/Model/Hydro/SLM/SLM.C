//******************************************************************************
/*!
  \file      src/Model/Hydro/SLM/SLM.C
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 06:12:45 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************

#include <Hydro/SLM/SLM.h>

using quinoa::SLM;

// void
// SLM::init()
// //******************************************************************************
// //  Initialize particle properties
// //! \author  J. Bakosi
// //******************************************************************************
// {
// //  initGaussian(m_velocities, m_nvelocity, rndstr(), m_str[0], 0.0, 1.0);
// //   tk::real r[m_nvelocity];
// // 
// //   // Generate initial values for all scalars for all particles
// //   for (uint64_t p=0; p<m_npar; ++p) {
// //     m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
// //                        m_str[0], m_nvelocity, r, 0.0, 1.0);
// //     memcpy(m_velocities + p*m_nvelocity, r, m_nvelocity*sizeof(tk::real));
// //   }
// }
// 
// void
// SLM::advance(int p, int tid, tk::real dt)
// //******************************************************************************
// //  Advance particles with the simplified Langevin model
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
