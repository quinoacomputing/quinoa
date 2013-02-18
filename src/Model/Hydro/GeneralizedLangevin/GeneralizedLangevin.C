//******************************************************************************
/*!
  \file      src/Model/Hydro/GeneralizedLangevin/GeneralizedLangevin.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 10:11:15 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydrodynamics model
  \details   Generalized Langevin hydrodynamics model
*/
//******************************************************************************

#include <iostream>

#include <GeneralizedLangevin.h>
#include <Hydro.h>

using namespace std;
using namespace Quinoa;

GeneralizedLangevin::GeneralizedLangevin(Memory* memory,
                                         Paradigm* paradigm) :
  Hydro(memory, paradigm, "Generalized Langevin")
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \author  J. Bakosi
//******************************************************************************
{
}

void
GeneralizedLangevin::echo()
//******************************************************************************
//  Echo information on the generalized Langevin model
//! \author  J. Bakosi
//******************************************************************************
{
}
