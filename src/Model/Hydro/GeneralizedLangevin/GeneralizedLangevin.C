//******************************************************************************
/*!
  \file      src/Model/Hydro/GeneralizedLangevin/GeneralizedLangevin.C
  \author    J. Bakosi
  \date      Fri Apr 26 16:36:21 2013
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
using namespace Quinoa;

GeneralizedLangevin::GeneralizedLangevin(Memory* const memory,
                                         Paradigm* const paradigm,
                                         Control* const control) :
  Hydro(memory, paradigm, control, "Generalized Langevin")
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \author  J. Bakosi
//******************************************************************************
{
}

void
GeneralizedLangevin::echo() const
//******************************************************************************
//  Echo information on the generalized Langevin model
//! \author  J. Bakosi
//******************************************************************************
{
}

void
GeneralizedLangevin::init()
//******************************************************************************
//  Initialize scalars
//! \author  J. Bakosi
//******************************************************************************
{
}

void
GeneralizedLangevin::advance(const real dt)
//******************************************************************************
//  Advance particles with the generalized Langevin model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  IGNORE(dt);
}
