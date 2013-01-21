//******************************************************************************
/*!
  \file      src/Model/HydroModel/HydroModel.C
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 10:31:53 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model base
  \details   Hydro model base
*/
//******************************************************************************

#include <HydroModel.h>
#include <HydroModelException.h>
#include <MemoryException.h>

using namespace Quinoa;

HydroModel::HydroModel(const string& name) :
  m_name(name)
//******************************************************************************
//  Constructor
//! \param[in]  name     Hydro model name
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_name.size() > 0, MemoryException,FATAL,EMPTY_NAME);
}
