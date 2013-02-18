//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 01:27:19 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model base
  \details   Hydro model base
*/
//******************************************************************************

#include <Hydro.h>
#include <HydroException.h>

using namespace Quinoa;

Hydro::Hydro(Memory* const memory,
             Paradigm* const paradigm,
             Control* const control,
             const string& name) :
  Model(memory, paradigm, control, name)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  name     Hydro model name
//! \author  J. Bakosi
//******************************************************************************
{
}
