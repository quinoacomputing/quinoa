//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.C
  \author    J. Bakosi
  \date      Sun 30 Sep 2012 09:23:35 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Unstructured mesh class definition
  \details   Unstructured mesh class definition
*/
//******************************************************************************

#include <iostream>

#include <UnsMesh.h>

using namespace Quinoa;

UnsMesh::~UnsMesh()
//******************************************************************************
//  Destructor: graceful free memory entries held
//! \author J. Bakosi
//******************************************************************************
{
  if (m_connLine) m_memory->freeEntry(m_connLine);
  if (m_connTri)  m_memory->freeEntry(m_connTri);
}
