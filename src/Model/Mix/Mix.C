//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 01:41:52 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

#include <Mix.h>
#include <MixException.h>
#include <Control.h>

using namespace Quinoa;
using namespace control;

Mix::Mix(Memory* const memory,
         Paradigm* const paradigm,
         Control* const control,
         const string& name) :
  Model(memory, paradigm, control, name),
  m_nscalar(control->get<NSCALAR>()),
  m_npar(control->get<NPAR>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  name     Mix model name
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_nscalar > 0, MixException,FATAL,BAD_NSCALAR);
  Assert(m_npar > 0, MixException,FATAL,BAD_NPAR);
}
