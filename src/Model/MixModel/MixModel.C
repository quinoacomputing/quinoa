//******************************************************************************
/*!
  \file      src/Model/MixModel/MixModel.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 01:31:57 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

#include <MixModel.h>
#include <MixModelException.h>

using namespace Quinoa;

MixModel::MixModel(Model* model, const string& name, const int& nscalar) :
  m_model(model), m_name(name), m_nscalar(nscalar)
//******************************************************************************
//  Constructor
//! \param[in]  model    Model object pointer
//! \param[in]  name     Mix model name
//! \param[in]  nscalar  Number of mixing scalars
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_nscalar > 0, MixModelException,FATAL,BAD_NSCALARS);
}
