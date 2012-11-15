//******************************************************************************
/*!
  \file      src/Model/MixModel/MixModel.C
  \author    J. Bakosi
  \date      Thu Nov 15 15:44:18 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

#include <MixModel.h>
#include <MixModelException.h>
#include <MemoryException.h>

using namespace Quinoa;

MixModel::MixModel(const int& nscalar, const string& name) :
  m_nscalar(nscalar), m_name(name)
//******************************************************************************
//  Constructor
//! \param[in]  nscalar  Number of mixing scalars
//! \param[in]  name     Mix model name
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_nscalar > 0, MixModelException,FATAL,BAD_NSCALARS);
  Assert(m_name.size() > 0, MemoryException,FATAL,EMPTY_NAME);
}
