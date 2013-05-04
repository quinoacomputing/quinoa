//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.C
  \author    J. Bakosi
  \date      Sat 04 May 2013 06:31:41 AM MDT
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
  if (m_nscalar <= 0)
    throw Exception(FATAL, "Wrong number of scalars: ", m_nscalar);

  if (m_npar <= 0)
    throw Exception(FATAL, "Wrong number of particles: ", m_npar);
}
