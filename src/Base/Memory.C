//******************************************************************************
/*!
  \file      src/Base/Memory.C
  \author    J. Bakosi
  \date      Sun 02 Sep 2012 03:00:26 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory (a store for MemoryEntry objects) base class declaration
  \details   Memory (a store for MemoryEntry objects) base class declaration
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
  m_entries = 0;
}

Memory::~Memory()
//******************************************************************************
//  Destructor
//! \details Free all allocated memory when leaving scope
//! \author  J. Bakosi
//******************************************************************************
{
  freeAll();
}

Int
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
  m_entry.push_back(entry);

  // Increase total number of memory entries, return id to caller
  return m_entries++;
}

void
Memory::freeAll()
//******************************************************************************
//  Deallocate all memory entries
//! \author J. Bakosi
//******************************************************************************
{
  if (m_entries) {
    do delete m_entry[--m_entries]; while (m_entries);
    m_entry.clear();
  }
}

size_t
Memory::getBytes()
//******************************************************************************
//  Get number of allocated bytes
//! \return Number of allocated bytes in the memory store
//! \author J. Bakosi
//******************************************************************************
{
  typedef vector<MemoryEntry*>::size_type size_type;

  if (m_entries) {
    size_t bytes = 0;
    size_type size = m_entry.size();
    for (size_type i=0; i<size; i++) {
      MemoryEntry* e = m_entry[i];
      bytes += sizeof(MemoryEntry) +
           e->m_number * VariableComponents[e->m_variable] * SizeOf[e->m_value];
    }
    return bytes;
  } else {
    return 0;
  }
}
