//******************************************************************************
/*!
  \file      src/Mesh/Mesh.C
  \author    J. Bakosi
  \date      Tue 11 Sep 2012 04:56:44 PM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh base class definition
  \details   Mesh base class definition
*/
//******************************************************************************

#include <algorithm>
#include <iterator>
#include <iostream>

#include <Mesh.h>
#include <Memory.h>
#include <MemoryException.h>

using namespace Quinoa;

Mesh::~Mesh()
//******************************************************************************
//  Destructor: free sets
//! \author J. Bakosi
//******************************************************************************
{
  // Free node sets
  if (m_entry.size()) {
    // Free all node sets
    MeshSet::const_iterator it;
    for (it=m_entry.begin(); it!=m_entry.end(); it++) {
      m_memory->freeEntry(*it);
    }
    // Clear container
    m_entry.clear();
  }
}

// void
// Mesh::newEntry(size_t number,
//                ValueType value,
//                VariableType variable,
//                string name,
//                Bool plot = false,
//                Bool restart = false)
// //******************************************************************************
// //  Add new MeshSet entry (e.g. list of nodes, elements, node ids, etc.)
// //! \author J. Bakosi
// //******************************************************************************
// {
//   MemoryEntry* entry = m_memory->newEntry(number,
//                                           value,
//                                           variable,
//                                           name,
//                                           plot,
//                                           restart);
// 
//   // Store new entry
//   pair<MeshSet::iterator,Bool> n = m_entry.insert(entry);
//   if (!n.second) throw MemoryException(FATAL, BAD_INSERT);
// }
