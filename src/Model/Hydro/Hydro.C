//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.C
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 04:30:48 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model base
  \details   Hydro model base
*/
//******************************************************************************

#include <Hydro.h>
#include <HydroException.h>
#include <Control.h>

using namespace Quinoa;

Hydro::Hydro(Memory* const memory,
             Paradigm* const paradigm,
             Control* const control,
             const string& name) :
  Model(memory, paradigm, control, name),
  m_nprop(NPROP),
  m_npar(control->get<control::NPAR>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  name     Hydro model name
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_nprop > 0, HydroException,FATAL,BAD_NPROP);
  Assert(m_npar > 0, ModelException,FATAL,BAD_NPAR);
}
