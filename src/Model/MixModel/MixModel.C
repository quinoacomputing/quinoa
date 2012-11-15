//******************************************************************************
/*!
  \file      src/Model/MixModel/MixModel.C
  \author    J. Bakosi
  \date      Thu Nov 15 13:33:23 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model base
  \details   Mix model base
*/
//******************************************************************************

#include <MixModel.h>
#include <MixModelException.h>

using namespace Quinoa;

MixModel::MixModel(Model* model,
                   MKLRandom* random,
                   const string& name,
                   const int& nscalar) :
  m_model(model), m_random(random), m_name(name), m_nscalar(nscalar)
//******************************************************************************
//  Constructor
//! \param[in]  model    Model object pointer
//! \param[in]  random   Random number generator object pointer
//! \param[in]  name     Mix model name
//! \param[in]  nscalar  Number of mixing scalars
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_nscalar > 0, MixModelException,FATAL,BAD_NSCALARS);
}
