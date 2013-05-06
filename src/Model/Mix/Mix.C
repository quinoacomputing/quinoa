//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.C
  \author    J. Bakosi
  \date      Mon May  6 11:07:47 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

#include <Mix.h>
#include <Control.h>
#include <Exception.h>

using namespace Quinoa;

Mix::Mix(Memory* const memory,
         Paradigm* const paradigm,
         Control* const control,
         const string& name) :
  Model(memory, paradigm, control, name),
  m_nscalar(control->get<control::NSCALAR>()),
  m_npar(control->get<control::NPAR>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  name     Mix model name
//! \author  J. Bakosi
//******************************************************************************
{
  Errchk(m_nscalar > 0, FATAL, "Wrong number of scalars");
  Errchk(m_npar > 0, FATAL, "Wrong number of particles");
}
