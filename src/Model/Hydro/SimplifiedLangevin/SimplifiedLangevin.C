//******************************************************************************
/*!
  \file      src/Model/Hydro/SimplifiedLangevin/SimplifiedLangevin.C
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 12:01:02 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************

#include <iostream>

#include <SimplifiedLangevin.h>
#include <Hydro.h>

using namespace std;
using namespace Quinoa;

SimplifiedLangevin::SimplifiedLangevin(Memory* const memory,
                                       Paradigm* const paradigm,
                                       Control* const control) :
  Hydro(memory, paradigm, control, "Simplified Langevin")
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
SimplifiedLangevin::echo()
//******************************************************************************
//  Echo information on the simplified Langevin model
//! \author  J. Bakosi
//******************************************************************************
{
}

void
SimplifiedLangevin::init()
//******************************************************************************
//  Initialize scalars
//! \author  J. Bakosi
//******************************************************************************
{
}

void
SimplifiedLangevin::advance(const real dt)
//******************************************************************************
//  Advance particles with the simplified Langevin model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
}
