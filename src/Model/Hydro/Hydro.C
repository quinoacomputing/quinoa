//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.C
  \author    J. Bakosi
  \date      Mon May  6 11:09:33 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model base
  \details   Hydro model base
*/
//******************************************************************************

#include <Hydro.h>
#include <Control.h>
#include <Exception.h>

using namespace Quinoa;

Hydro::Hydro(Memory* const memory,
             Paradigm* const paradigm,
             Control* const control,
             const string& name) :
  Model(memory, paradigm, control, name),
  m_nprop(6),
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
  Errchk(m_nprop > 0, FATAL, "Wrong number of particle properties");
  Errchk(m_npar > 0, FATAL, "Wrong number of particles");
}
