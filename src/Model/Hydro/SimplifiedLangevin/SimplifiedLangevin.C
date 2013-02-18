//******************************************************************************
/*!
  \file      src/Model/Hydro/SimplifiedLangevin/SimplifiedLangevin.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 10:10:08 AM MST
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

SimplifiedLangevin::SimplifiedLangevin(Memory* memory,
                                       Paradigm* paradigm) :
  Hydro(memory, paradigm, "Simplified Langevin")
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
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
