//******************************************************************************
/*!
  \file      src/Model/Hydro/SLM/SLM.C
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 08:37:34 PM MDT
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
#include <Hydro/Hydro.h>
#include <Hydro/SLM/SLM.h>

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
SimplifiedLangevin::advance(int p, int tid, real dt)
//******************************************************************************
//  Advance particles with the simplified Langevin model
//! \param[in]  p    Particle to advance
//! \param[in]  tid  Thread id
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
IGNORE(p);
IGNORE(tid);
IGNORE(dt);
}
