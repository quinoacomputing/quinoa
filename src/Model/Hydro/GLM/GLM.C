//******************************************************************************
/*!
  \file      src/Model/Hydro/GLM/GLM.C
  \author    J. Bakosi
  \date      Fri Sep 27 11:59:38 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydrodynamics model
  \details   Generalized Langevin hydrodynamics model
*/
//******************************************************************************

#include <iostream>

#include <Macro.h>
#include <GLM.h>
#include <Hydro.h>

using namespace std;
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
