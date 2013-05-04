//******************************************************************************
/*!
  \file      src/Model/Model.C
  \author    J. Bakosi
  \date      Fri 03 May 2013 06:17:13 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************

#include <iostream>
#include <cassert>

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
  assert(m_name.size() > 0);
}
