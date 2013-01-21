//******************************************************************************
/*!
  \file      src/Model/Mix/Mix.C
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 11:52:44 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

#include <Mix.h>
#include <MixException.h>

using namespace Quinoa;

Mix::Mix(const int& nscalar, const string& name) :
  Model(name),
  m_nscalar(nscalar)
//******************************************************************************
//  Constructor
//! \param[in]  nscalar  Number of mixing scalars
//! \param[in]  name     Mix model name
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_nscalar > 0, MixException,FATAL,BAD_NSCALARS);
}
