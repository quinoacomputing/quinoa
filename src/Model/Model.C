//******************************************************************************
/*!
  \file      src/Model/Model.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 01:17:14 PM MST
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

Model::Model(Memory* const memory,
             Paradigm* const paradigm,
             Control* const control,
             const string& name) :
  m_memory(memory), m_paradigm(paradigm), m_control(control), m_name(name)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  name     Name of model
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_name.size() > 0, MemoryException,FATAL,EMPTY_NAME);
}
