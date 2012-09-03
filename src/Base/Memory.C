//******************************************************************************
/*!
  \file      src/Base/Memory.C
  \author    J. Bakosi
  \date      Sun 02 Sep 2012 07:00:11 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory (a store for MemoryEntry objects) base class definition
  \details   Memory (a store for MemoryEntry objects) base class definition
*/
//******************************************************************************

#include <string>
#include <cassert>
#include <iostream>

using namespace std;

#include <MemoryEntry.h>
#include <Memory.h>

using namespace Quinoa;

Memory::Memory()
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
}

Memory::~Memory()
//******************************************************************************
//  Destructor
//! \details Free all allocated memory when leaving scope (just in case)
//! \author  J. Bakosi
//******************************************************************************
{
  freeAllEntries();
}

MemoryEntry*
Memory::newEntry(size_t number,
                 ValueType value,
                 VariableType variable,
                 string name,
                 Bool plot,
                 Bool restart)
//******************************************************************************
//  Allocate memory entry
//! \param[in]  number    Number of items to allocate
//! \param[in]  value     Chosen value type (BOOL_VAL, INT_VAL, REAL_VAL, etc.)
//! \param[in]  variable  Chosen variable type (SCALAR_VAR, VECTOR_VAR, etc.)
//! \param[in]  name      Given variable name
//! \param[in]  plot      True if variable can be plotted
//! \param[in]  restart   True if variable will be written to restart file
//! \return               ID for the allocated entry
//! \author     J. Bakosi
//******************************************************************************
{
  assert(number > 0);
  assert(0 <= value && value < NUM_VALUE_TYPES);
  assert(0 <= variable && variable < NUM_VARIABLE_TYPES);
  assert(name.size() > 0);

  // Compute total number of bytes to be allocated
  size_t nbytes = number * VariableComponents[variable] * SizeOf[value];

  // Allocate memory
  void* ptr = static_cast<void*>(new char [nbytes]);

  // Allocate memory entry
  MemoryEntry* entry = new MemoryEntry(number,
                                       value,
                                       variable,
                                       name,
                                       plot,
                                       restart,
                                       ptr);
  // Store memory entry
  pair<MemorySet::iterator,Bool> p = m_entry.insert(entry);
  cout << p.second << ": " << ptr << endl;

  // Return key to caller
  return entry;
}

void
Memory::freeEntry(MemoryEntry* id)
//******************************************************************************
//  Deallocate a memory entry
//! \param[in]  id  ID of the entry to be freed
//! \author J. Bakosi
//******************************************************************************
{
  // Return if entry is already deallocated
  if (id == 0) return;

  if (m_entry.size()) {
    // Deallocate memory entry pointed to by m_entry[id]
    // This also automatically calls MemoryEntry::~MemoryEntry(), which
    // deallocates the memory pointed to by MemoryEntry::m_ptr
    auto it = m_entry.find(id);
    if (it!=m_entry.end()) {
      delete *it;
    }

    // Remove MemoryEntry from MemorySet
    m_entry.erase(id);

    // Zero id, so the caller can also tell that this memory entry has been
    // deallocated
    id = 0;
  }
}

void
Memory::freeAllEntries()
//******************************************************************************
//  Deallocate all memory entries
//! \author J. Bakosi
//******************************************************************************
{
  if (m_entry.size()) {
    for (auto it=m_entry.begin(); it!=m_entry.end(); it++) {
      delete *it;
    }
    m_entry.clear();
  }
}

size_t
Memory::getBytes()
//******************************************************************************
//  Get number of allocated bytes
//! \details Return the number of bytes allocated in newEntry(). We account for
//!          the size of the MemoryEntry class instances and the allocated data
//!          pointed to by MemoryEntry::m_ptr. We do not account for the
//!          overhead of the STL container, therefore we will always
//!          underestimate the actual memory usage, though by only a very small
//!          fraction, i.e. <1e-4% for memory allocated in the range of MBytes.
//! \return Number of allocated bytes in the memory store
//! \author J. Bakosi
//******************************************************************************
{
  size_t bytes = 0;

  if (m_entry.size()) {
    for (auto it=m_entry.begin(); it!=m_entry.end(); it++) {
      bytes += sizeof(MemoryEntry) +
                 (*it)->m_number *
                 VariableComponents[(*it)->m_variable] *
                 SizeOf[(*it)->m_value];
     }
  }

  return bytes;
}
