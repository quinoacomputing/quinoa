//******************************************************************************
/*!
  \file      src/Model/Model.C
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 11:53:39 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************

#include <iostream>

#include <MemoryException.h>
#include <Model.h>

using namespace std;
using namespace Quinoa;

Model::Model(const string& name) : m_name(name)
//******************************************************************************
//  Constructor
//! \param[in]  name     Name of model
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_name.size() > 0, MemoryException,FATAL,EMPTY_NAME);
}

Model::~Model()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
}
