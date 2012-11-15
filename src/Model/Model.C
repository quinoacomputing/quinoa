//******************************************************************************
/*!
  \file      src/Model/Model.C
  \author    J. Bakosi
  \date      Thu Nov 15 15:17:52 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************

#include <iostream>

#include <Memory.h>
#include <Model.h>

using namespace std;
using namespace Quinoa;

Model::Model(Memory* memory, Paradigm* paradigm, const string& name) :
  m_memory(memory), m_paradigm(paradigm), m_name(name)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  name     Name of model
//! \author  J. Bakosi
//******************************************************************************
{
}

Model::~Model()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
}
