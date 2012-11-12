//******************************************************************************
/*!
  \file      src/Model/MixModel/MixModel.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 09:12:19 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

#include <cassert>

#include <MixModel.h>

using namespace Quinoa;

MixModel::MixModel(const int& nscalar) : m_nscalar(nscalar)
//******************************************************************************
//  Constructor
//! \param[in]  nscalar  Number of mixing scalars
//! \author  J. Bakosi
//******************************************************************************
{
  assert(m_nscalar > 0);
}
