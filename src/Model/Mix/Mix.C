//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 09:57:14 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

#include <Mix.h>
#include <MixException.h>

using namespace Quinoa;

Mix::Mix(Memory* memory,
         Paradigm* paradigm,
         const int& nscalar,
         const int& npar,
         const string& name) :
  Model(memory, paradigm, name),
  m_nscalar(nscalar), m_npar(npar)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  npar     Number of particles
//! \param[in]  nscalar  Number of mixing scalars
//! \param[in]  name     Mix model name
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_nscalar > 0, MixException,FATAL,BAD_NSCALAR);
  Assert(m_npar > 0, MixException,FATAL,BAD_NPAR);
}
