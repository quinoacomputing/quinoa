//******************************************************************************
/*!
  \file      src/Model/Hydro/GLM/GLM.C
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 08:37:26 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydrodynamics model
  \details   Generalized Langevin hydrodynamics model
*/
//******************************************************************************

#include <iostream>

#include <Macro.h>
#include <Hydro/Hydro.h>
#include <Hydro/GLM/GLM.h>

using namespace quinoa;

void
GeneralizedLangevin::init()
//******************************************************************************
//  Initialize scalars
//! \author  J. Bakosi
//******************************************************************************
{
}

void
GeneralizedLangevin::advance(int p, int tid, real dt)
//******************************************************************************
//  Advance particles with the generalized Langevin model
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
