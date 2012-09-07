//******************************************************************************
/*!
  \file      src/Base/Memory.C
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 03:35:12 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory (a store for MemoryEntry objects) base class definition
  \details   Memory (a store for MemoryEntry objects) base class definition
*/
//******************************************************************************

#include <cassert>
#include <cstring>

using namespace std;

#include <MemoryEntry.h>
#include <Memory.h>
#include <MemoryException.h>

using namespace Quinoa;

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

  // Allocate memory and its metadata
  try {
    void* ptr = static_cast<void*>(new char [nbytes]);

    MemoryEntry* entry = new MemoryEntry(nbytes,
                                         number,
                                         value,
                                         variable,
                                         name,
                                         plot,
                                         restart,
                                         ptr);
    // Store memory entry
    pair<MemorySet::iterator,Bool> e = m_entry.insert(entry);
    if (!e.second) throw MemoryException(FATAL, BAD_INSERT);

    // Map variable name to MemorySet key
    pair<MemoryNames::iterator,Bool> n = m_name.emplace(name,entry);
    if (!n.second) throw MemoryException(FATAL, BAD_INSERT);

    // Test if name is unique
    if (m_entry.size() != m_name.size())
      throw MemoryException(WARNING, BAD_NAME);

    // Return key to caller
    return entry;
  } catch (bad_alloc& ba) { throw MemoryException(FATAL, BAD_ALLOC); }
}

MemoryEntry*
Memory::newZeroEntry(size_t number,
                     ValueType value,
                     VariableType variable,
                     string name,
                     Bool plot,
                     Bool restart)
//******************************************************************************
//  Allocate and zero memory entry
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
  // Allocate new entry
  MemoryEntry* entry = newEntry(number, value, variable, name, plot, restart);

  // Zero new entry
  zero(entry);

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
  // Return and throw warning if entry is already deallocated
  if (id == 0) throw MemoryException(WARNING, UNDEFINED);

  // Return and throw warning if memory store is empty
  if (!m_entry.size()) throw MemoryException(WARNING, EMPTY_STORE);

  // Find and remove memory entry
  auto it = m_entry.find(id);
  if (it==m_entry.end()) throw MemoryException(WARNING, NOT_FOUND);

  // Remove variable name mapped to MemorySet key
  MemoryNames::size_type erased = m_name.erase((*it)->m_name);
  if (!erased) throw MemoryException(WARNING, NOT_ERASED);

  // Deallocate memory entry pointed to by m_entry[id]
  // This also automatically calls MemoryEntry::~MemoryEntry(), which
  // deallocates the memory pointed to by MemoryEntry::m_ptr
  delete *it;

  // Remove MemoryEntry from MemorySet
  MemorySet::size_type removed = m_entry.erase(id);
  if (!removed) throw MemoryException(WARNING, NOT_ERASED);

  // Zero id, so the caller can also tell that this memory entry has been
  // deallocated
  id = const_cast<MemoryEntry*>(UNDEFINED_ENTRY);
}

void
Memory::freeAllEntries() noexcept
//******************************************************************************
//  Deallocate all memory entries
//! \details This may (and probably will) be called several times. At normal
//!          leave-of-scope the destructor calls it. The Driver also calls it
//!          via its destructor. Additionally, the Exception::handleException()
//!          may also call it (through Driver::finalize()) at any time if a
//!          FATAL error is encountered. Because of this the test of
//!          m_entry.size() does not throw an exception.
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
  // Return and throw warning if memory store is empty
  if (!m_entry.size()) throw MemoryException(WARNING, EMPTY_STORE);

  // Find memory entry and return its number of variables
  auto it = m_entry.find(id);
  if (it==m_entry.end()) throw MemoryException(WARNING, NOT_FOUND);
  return (*it)->m_number;
}

ValueType
Memory::getValue(MemoryEntry* id)
//******************************************************************************
//  Return the value type based on the ID
//! \return The value type
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if memory store is empty
  if (!m_entry.size()) throw MemoryException(WARNING, EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  if (it==m_entry.end()) throw MemoryException(WARNING, NOT_FOUND);
  return (*it)->m_value;
}

VariableType
Memory::getVariable(MemoryEntry* id)
//******************************************************************************
//  Return the variable type based on the ID
//! \return The variable type
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if memory store is empty
  if (!m_entry.size()) throw MemoryException(WARNING, EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  if (it==m_entry.end()) throw MemoryException(WARNING, NOT_FOUND);
  return (*it)->m_variable;
}

string
Memory::getName(MemoryEntry* id)
//******************************************************************************
//  Return the variable name based on the ID
//! \return The variable name
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if memory store is empty
  if (!m_entry.size()) throw MemoryException(WARNING, EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  if (it==m_entry.end()) throw MemoryException(WARNING, NOT_FOUND);
  return (*it)->m_name;
}

Bool
Memory::getPlot(MemoryEntry* id)
//******************************************************************************
//  Return true if the variable can be plotted based on the ID
//! \return True if the variable can be plotted
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if memory store is empty
  if (!m_entry.size()) throw MemoryException(WARNING, EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  if (it==m_entry.end()) throw MemoryException(WARNING, NOT_FOUND);
  return (*it)->m_plot;
}

Bool
Memory::getRestart(MemoryEntry* id)
//******************************************************************************
//  Return true if the variable is writted to restart file based on the ID
//! \return True if the variable is written to restart file
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if memory store is empty
  if (!m_entry.size()) throw MemoryException(WARNING, EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  if (it==m_entry.end()) throw MemoryException(WARNING, NOT_FOUND);
  return (*it)->m_restart;
}

MemoryEntry*
Memory::getID(string name)
//******************************************************************************
//  Return the MemorySet key based on the variable name
//! \return MemorySet key (used as ID)
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if memory store is empty
  if (!m_name.size()) throw MemoryException(WARNING, EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_name.find(name);
  if (it==m_name.end()) throw MemoryException(WARNING, NOT_FOUND);
  return it->second;
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
  // Return and throw warning if memory store is empty
  if (!m_entry.size()) throw MemoryException(WARNING, EMPTY_STORE);

  size_t bytes = 0;
  for (auto it=m_entry.begin(); it!=m_entry.end(); it++) {
    bytes += sizeof(MemoryEntry) +
               (*it)->m_number * VariableComponents[(*it)->m_variable] *
               SizeOf[(*it)->m_value];
  }
  return bytes;
}

void
Memory::zero(MemoryEntry* id)
//******************************************************************************
//  Zero entry using multiple threads
//! \author J. Bakosi
//******************************************************************************
{
  // Get size of value type
  size_t size = SizeOf[id->m_value];

  // Compute chunk size
  Int i = id->m_number/m_nthreads;

  // Zero remaining portion
  memset(static_cast<char*>(id->m_ptr) + m_nthreads*i*size,
         0,
         (id->m_number%m_nthreads)*size);

  Int myid;
  #ifdef _OPENMP
  #pragma omp parallel private(myid)
  #endif
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif

    // each processor zeros its own portion
    memset(static_cast<char*>(id->m_ptr) + myid*i*size, 0, i*size);
  }
}
