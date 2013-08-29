//******************************************************************************
/*!
  \file      src/Model/Hydro/GeneralizedLangevin/GeneralizedLangevin.C
  \author    J. Bakosi
  \date      Thu Aug 29 15:34:46 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydrodynamics model
  \details   Generalized Langevin hydrodynamics model
*/
//******************************************************************************

#include <iostream>

#include <Macro.h>
#include <GeneralizedLangevin.h>
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
GeneralizedLangevin::advance(const real& dt)
//******************************************************************************
//  Advance particles with the generalized Langevin model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  IGNORE(dt);
}
