//******************************************************************************
/*!
  \file      src/Model/Hydro/GeneralizedLangevin/GeneralizedLangevin.C
  \author    J. Bakosi
  \date      Mon 13 May 2013 09:48:10 PM MDT
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
                                         Control* const control,
                                         real* const velocities) :
  Hydro<GeneralizedLangevin>(memory,
                             paradigm,
                             control,
                             control->get<control::NVELOCITY>(),
                             velocities)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  velocities Pointer to particle velocities
//! \author  J. Bakosi
//******************************************************************************
{
  IGNORE(m_velocities);
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
GeneralizedLangevin::advance(const real& dt)
//******************************************************************************
//  Advance particles with the generalized Langevin model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  IGNORE(dt);
}
