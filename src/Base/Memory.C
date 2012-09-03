//******************************************************************************
/*!
  \file      src/Base/Memory.C
  \author    J. Bakosi
  \date      Sun 02 Sep 2012 11:52:52 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory (a store for MemoryEntry objects) base class definition
  \details   Memory (a store for MemoryEntry objects) base class definition
*/
//******************************************************************************

#include <iostream>
#include <cassert>

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

  // Map variable name to MemorySet key
  m_name.emplace(name,entry);

  if (m_entry.size() != m_name.size()) {
    cout << "Names should be unique!" << endl;
  }

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
    auto it = m_entry.find(id);
    if (it!=m_entry.end()) {
      // Remove variable name mapped to MemorySet key
      m_name.erase((*it)->m_name);
      // Deallocate memory entry pointed to by m_entry[id]
      // This also automatically calls MemoryEntry::~MemoryEntry(), which
      // deallocates the memory pointed to by MemoryEntry::m_ptr
      delete *it;
    }

    // Remove MemoryEntry from MemorySet
    m_entry.erase(id);

    // Zero id, so the caller can also tell that this memory entry has been
    // deallocated
    id = const_cast<MemoryEntry*>(UNDEFINED_ENTRY);
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
    // Deallocate memory entry pointed to by m_entry[*]
    // This also automatically calls MemoryEntry::~MemoryEntry(), which
    // deallocates the memory pointed to by MemoryEntry::m_ptr
    for (auto it=m_entry.begin(); it!=m_entry.end(); it++) {
      delete *it;
    }
    // Clear containers
    m_name.clear();
    m_entry.clear();
  }
}

size_t
Memory::getNumber(MemoryEntry* id)
//******************************************************************************
//  Return the number of items based on the ID
//! \return Number of items
//! \author J. Bakosi
//******************************************************************************
{
  size_t number = 0;

  if (m_entry.size()) {
    auto it = m_entry.find(id);
    if (it!=m_entry.end())
      number = (*it)->m_number;
  }

  return number;
}

ValueType
Memory::getValue(MemoryEntry* id)
//******************************************************************************
//  Return the value type based on the ID
//! \return The value type
//! \author J. Bakosi
//******************************************************************************
{
  ValueType value;

  if (m_entry.size()) {
    auto it = m_entry.find(id);
    if (it!=m_entry.end())
      value = (*it)->m_value;
  }

  return value;
}

VariableType
Memory::getVariable(MemoryEntry* id)
//******************************************************************************
//  Return the variable type based on the ID
//! \return The variable type
//! \author J. Bakosi
//******************************************************************************
{
  VariableType variable;

  if (m_entry.size()) {
    auto it = m_entry.find(id);
    if (it!=m_entry.end())
      variable = (*it)->m_variable;
  }

  return variable;
}

string
Memory::getName(MemoryEntry* id)
//******************************************************************************
//  Return the variable name based on the ID
//! \return The variable name
//! \author J. Bakosi
//******************************************************************************
{
  string name;

  if (m_entry.size()) {
    auto it = m_entry.find(id);
    if (it!=m_entry.end())
      name = (*it)->m_name;
  }

  return name;
}

Bool
Memory::getPlot(MemoryEntry* id)
//******************************************************************************
//  Return true if the variable can be plotted based on the ID
//! \return True if the variable can be plotted
//! \author J. Bakosi
//******************************************************************************
{
  Bool plot;

  if (m_entry.size()) {
    auto it = m_entry.find(id);
    if (it!=m_entry.end())
      plot = (*it)->m_plot;
  }

  return plot;
}

Bool
Memory::getRestart(MemoryEntry* id)
//******************************************************************************
//  Return true if the variable is writted to restart file based on the ID
//! \return True if the variable is written to restart file
//! \author J. Bakosi
//******************************************************************************
{
  Bool restart;

  if (m_entry.size()) {
    auto it = m_entry.find(id);
    if (it!=m_entry.end())
      restart = (*it)->m_restart;
  }

  return restart;
}

MemoryEntry*
Memory::getID(string name)
//******************************************************************************
//  Return the MemorySet key based on the variable name
//! \return MemorySet key (used as ID)
//! \author J. Bakosi
//******************************************************************************
{
  auto it = m_name.find(name);
  if (it!=m_name.end())
    return it->second;
  else
    return 0;
}

size_t
Memory::getBytes()
//******************************************************************************
//  Return the number of allocated bytes
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
