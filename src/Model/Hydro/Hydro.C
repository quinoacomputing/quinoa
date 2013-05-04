//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.C
  \author    J. Bakosi
  \date      Sat 04 May 2013 06:33:50 AM MDT
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
  if (m_nprop <= 0)
    throw Exception(FATAL, "Wrong number of particle properties: ", m_nprop);

  if (m_npar <= 0)
    throw Exception(FATAL, "Wrong number of particles: ", m_npar);
}
