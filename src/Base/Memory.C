//******************************************************************************
/*!
  \file      src/Base/Memory.C
  \author    J. Bakosi
  \date      Sun 11 Nov 2012 12:43:48 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Memory (a store for MemoryEntry objects) base class definition
  \details   Memory (a store for MemoryEntry objects) base class definition
*/
//******************************************************************************

#include <cstring>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

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
                 ValType value,
                 VarType variable,
                 string name,
                 bool plot,
                 bool restart)
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
  Assert(number > 0, MemoryException,FATAL,BAD_NUMBER);
  Assert(name.size() > 0, MemoryException,FATAL,EMPTY_NAME);

  // Compute total number of bytes to be allocated
  size_t nbytes = number *
                  VarComp[static_cast<int>(variable)] *
                  SizeOf[static_cast<int>(value)];

  // Allocate memory
  void* ptr = static_cast<void*>(new (nothrow) char [nbytes]);
  Assert(ptr != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Allocate memory entry metadata
  MemoryEntry* entry = new (nothrow) MemoryEntry(nbytes,
                                                 number,
                                                 value,
                                                 variable,
                                                 name,
                                                 plot,
                                                 restart,
                                                 ptr);
  if (entry == nullptr) {
    if (ptr) {
      delete [] static_cast<char*>(ptr);
      ptr = nullptr;
    }
    Throw(MemoryException,FATAL,BAD_ALLOC);
  }

  // Store memory entry
  pair<MemorySet::iterator,bool> e = m_entry.insert(entry);
  if (!e.second) {
    if (entry) {
       delete entry;
       entry = nullptr;
    }
    Throw(MemoryException,FATAL,BAD_INSERT);
  }

  // Map variable name to MemorySet key
  pair<MemoryNames::iterator,bool> n = m_name.emplace(name,entry);
  if (!n.second) {
    if (entry) {
      delete entry;
      entry = nullptr;
    }
    Throw(MemoryException,FATAL,BAD_INSERT);
  }

  // Test if name is unique
  Assert(m_entry.size()==m_name.size(), MemoryException,WARNING,NONUNIQUE_NAME);

  // Return key to caller
  return entry;
}

MemoryEntry*
Memory::newZeroEntry(size_t number,
                     ValType value,
                     VarType variable,
                     string name,
                     bool plot,
                     bool restart)
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
  Assert(id != nullptr, MemoryException,WARNING,UNDEFINED);

  // Return and throw warning if memory store is empty
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Find memory entry
  auto it = m_entry.find(id);
  Assert(it != m_entry.end(), MemoryException,WARNING,NOT_FOUND);

  // Remove variable name mapped to MemorySet key
  Errchk(m_name.erase((*it)->m_name), MemoryException,WARNING,NOT_ERASED);

  // Deallocate memory entry pointed to by m_entry[id]
  // This also automatically calls MemoryEntry::~MemoryEntry(), which
  // deallocates the memory pointed to by MemoryEntry::m_ptr
  delete *it;

  // Remove MemoryEntry from MemorySet
  Errchk(m_entry.erase(id), MemoryException,WARNING,NOT_ERASED);

  // Zero id, so the caller can also tell that the memory entry has been removed
  id = nullptr;
}

void
Memory::freeAllEntries() noexcept
//******************************************************************************
//  Deallocate all memory entries
//! \details This may (and probably will) be called several times. At normal
//!          leave-of-scope the destructor calls it. The Driver also calls it
//!          via its destructor. Additionally, the Exception::handleException()
//!          may also call it (through Driver::finalize()) at any time if a
//!          FATAL error is encountered. Because of this we do not throw an
//!          exception here.
//! \author J. Bakosi
//******************************************************************************
{
  // Deallocate all memory entries
  // This also automatically calls MemoryEntry::~MemoryEntry(), which
  // deallocates the memory pointed to by MemoryEntry::m_ptr
  for (auto& e : m_entry) { delete e; }
  // Clear containers
  m_name.clear();
  m_entry.clear();
}

void
Memory::echoAllEntries(MemoryEntryField crit)
//******************************************************************************
//  Echo all memory entries
//! \param[in]  crit    Sort items by memory entry field criteria
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if memory store is empty
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Echo AllEntries-header
  cout << "* Dynamically allocated memory entries ";
  switch (crit) {
    case MemoryEntryField::BYTES:    cout << "(sorted by Bytes)\n";    break;
    case MemoryEntryField::NUMBER:   cout << "(sorted by Number)\n";   break;
    case MemoryEntryField::VALUE:    cout << "(sorted by Value)\n";    break;
    case MemoryEntryField::VARIABLE: cout << "(sorted by Variable)\n"; break;
    case MemoryEntryField::NAME:     cout << "(sorted by Name)\n";     break;
    case MemoryEntryField::PLOT:     cout << "(sorted by Plot)\n";     break;
    case MemoryEntryField::RESTART:  cout << "(sorted by Restart)\n";  break;
    case MemoryEntryField::UNSPECIFIED: default: cout << "\n";
  }
  cout << "  " << setw(EntryWidth[0]) << "Name"
       << "  " << setw(EntryWidth[1]) << "Number"
       << "  " << setw(EntryWidth[2]) << "Value"
       << "  " << setw(EntryWidth[3]) << "ValueSize"
       << "  " << setw(EntryWidth[4]) << "Variable"
       << "  " << setw(EntryWidth[5]) << "Bytes"
       << "  " << setw(EntryWidth[6]) << "Plot"
       << "  " << setw(EntryWidth[7]) << "Restart"
       << "  " << setw(EntryWidth[8]) << "Ptr"
       << "\n";
  cout << setfill('=');
  cout << "  " << setw(EntryWidth[0]) << "="
       << "  " << setw(EntryWidth[1]) << "="
       << "  " << setw(EntryWidth[2]) << "="
       << "  " << setw(EntryWidth[3]) << "="
       << "  " << setw(EntryWidth[4]) << "="
       << "  " << setw(EntryWidth[5]) << "="
       << "  " << setw(EntryWidth[6]) << "="
       << "  " << setw(EntryWidth[7]) << "="
       << "  " << setw(EntryWidth[8]) << "="
       << endl;
  cout << setfill(' ');

  switch (crit) {
    case MemoryEntryField::BYTES:    echoByBytes();    break;
    case MemoryEntryField::NUMBER:   echoByNumber();   break;
    case MemoryEntryField::VALUE:    echoByValue();    break;
    case MemoryEntryField::VARIABLE: echoByVariable(); break;
    case MemoryEntryField::NAME:     echoByName();     break;
    case MemoryEntryField::PLOT:     echoByPlot();     break;
    case MemoryEntryField::RESTART:  echoByRestart();  break;
    case MemoryEntryField::UNSPECIFIED: default: echo();
  }
}

void
Memory::echo()
//******************************************************************************
//  Echo unsorted memory entries
//! \author J. Bakosi
//******************************************************************************
{
  for (auto* e : m_entry) cout << e->line();
}

void
Memory::echoByBytes()
//******************************************************************************
//  Echo memory entries sorted by Bytes
//! \author J. Bakosi
//******************************************************************************
{
  // Copy unordered memory entry keys to vector
  vector<MemoryEntry*> srt;
  copy(m_entry.begin(), m_entry.end(), back_inserter(srt));

  // Sort vector of memory entries by their Bytes
  sort(srt.begin(), srt.end(),
       [] (const MemoryEntry* a, const MemoryEntry* b) {
         return a->m_bytes > b->m_bytes;
       });

  // Echo ordered entries
  for (auto* e : srt) cout << e->line();
}

void
Memory::echoByNumber()
//******************************************************************************
//  Echo memory entries sorted by Number
//! \author J. Bakosi
//******************************************************************************
{
  // Copy unordered memory entry keys to vector
  vector<MemoryEntry*> srt;
  copy(m_entry.begin(), m_entry.end(), back_inserter(srt));

  // Sort vector of memory entries by their Number
  sort(srt.begin(), srt.end(),
       [] (const MemoryEntry* a, const MemoryEntry* b) {
         return a->m_number > b->m_number;
       });

  // Echo ordered entries
  for (auto* e : srt) cout << e->line();
}


void
Memory::echoByValue()
//******************************************************************************
//  Echo memory entries sorted by Value
//! \author J. Bakosi
//******************************************************************************
{
  // Copy unordered memory entry keys to vector
  vector<MemoryEntry*> srt;
  copy(m_entry.begin(), m_entry.end(), back_inserter(srt));

  // Sort vector of memory entries by their Value
  sort(srt.begin(), srt.end(),
       [] (const MemoryEntry* a, const MemoryEntry* b) {
         return a->m_value < b->m_value;
       });

  // Echo ordered entries
  for (auto* e : srt) cout << e->line();
}

void
Memory::echoByVariable()
//******************************************************************************
//  Echo memory entries sorted by Variable
//! \author J. Bakosi
//******************************************************************************
{
  // Copy unordered memory entry keys to vector
  vector<MemoryEntry*> srt;
  copy(m_entry.begin(), m_entry.end(), back_inserter(srt));

  // Sort vector of memory entries by their Variable
  sort(srt.begin(), srt.end(),
       [] (const MemoryEntry* a, const MemoryEntry* b) {
         return a->m_variable < b->m_variable;
       });

  // Echo ordered entries
  for (auto* e : srt) cout << e->line();
}

void
Memory::echoByName()
//******************************************************************************
//  Echo memory entries sorted by Name
//! \author J. Bakosi
//******************************************************************************
{
  // Copy unordered memory entry keys to vector
  vector<MemoryEntry*> srt;
  copy(m_entry.begin(), m_entry.end(), back_inserter(srt));

  // Sort vector of memory entries by their Name
  sort(srt.begin(), srt.end(),
       [] (const MemoryEntry* a, const MemoryEntry* b) {
         return a->m_name < b->m_name;
       });

  // Echo ordered entries
  for (auto* e : srt) cout << e->line();
}

void
Memory::echoByPlot()
//******************************************************************************
//  Echo memory entries sorted by Plot
//! \author J. Bakosi
//******************************************************************************
{
  // Copy unordered memory entry keys to vector
  vector<MemoryEntry*> srt;
  copy(m_entry.begin(), m_entry.end(), back_inserter(srt));

  // Sort vector of memory entries by their Plot
  sort(srt.begin(), srt.end(),
       [] (const MemoryEntry* a, const MemoryEntry* b) {
         return a->m_plot < b->m_plot;
       });

  // Echo ordered entries
  for (auto* e : srt) cout << e->line();
}

void
Memory::echoByRestart()
//******************************************************************************
//  Echo memory entries sorted by Restart
//! \author J. Bakosi
//******************************************************************************
{
  // Copy unordered memory entry keys to vector
  vector<MemoryEntry*> srt;
  copy(m_entry.begin(), m_entry.end(), back_inserter(srt));

  // Sort vector of memory entries by their Restart
  sort(srt.begin(), srt.end(),
       [] (const MemoryEntry* a, const MemoryEntry* b) {
         return a->m_restart < b->m_restart;
       });

  // Echo ordered entries
  for (auto* e : srt) cout << e->line();
}

size_t
Memory::getNumber(MemoryEntry* id)
//******************************************************************************
//  Return the number of items based on the ID
//! \param[in]  id  ID of the entry
//! \return Number of items
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if entry is invalid
  Assert(id != nullptr, MemoryException,WARNING,UNDEFINED);

  // Return and throw warning if memory store is empty
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Find memory entry and return its number of variables
  auto it = m_entry.find(id);
  Assert(it!=m_entry.end(), MemoryException,WARNING,NOT_FOUND);
  return (*it)->m_number;
}

size_t
Memory::getNumber(const string& name)
//******************************************************************************
//  Return the number of items based on the variable name
//! \param[in]  name  name of the variable
//! \return Number of items
//! \author J. Bakosi
//******************************************************************************
{
  // Find memory entry and return its number of variables
  return getID(name)->m_number;
}

ValType
Memory::getValue(MemoryEntry* id)
//******************************************************************************
//  Return the value type based on the ID
//! \param[in]  id  ID of the entry
//! \return The value type
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if entry is invalid
  Assert(id != nullptr, MemoryException,WARNING,UNDEFINED);

  // Return and throw warning if memory store is empty
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  Assert(it != m_entry.end(), MemoryException,WARNING,NOT_FOUND);
  return (*it)->m_value;
}

VarType
Memory::getVariable(MemoryEntry* id)
//******************************************************************************
//  Return the variable type based on the ID
//! \param[in]  id  ID of the entry
//! \return The variable type
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if entry is invalid
  Assert(id != nullptr, MemoryException,WARNING,UNDEFINED);

  // Return and throw warning if memory store is empty
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  Assert(it != m_entry.end(), MemoryException,WARNING,NOT_FOUND);
  return (*it)->m_variable;
}

string
Memory::getName(MemoryEntry* id)
//******************************************************************************
//  Return the variable name based on the ID
//! \param[in]  id  ID of the entry
//! \return The variable name
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if entry is invalid
  Assert(id != nullptr, MemoryException,WARNING,UNDEFINED);

  // Return and throw warning if memory store is empty
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  Assert(it != m_entry.end(), MemoryException,WARNING,NOT_FOUND);
  return (*it)->m_name;
}

bool
Memory::getPlot(MemoryEntry* id)
//******************************************************************************
//  Return true if the variable can be plotted based on the ID
//! \param[in]  id  ID of the entry
//! \return True if the variable can be plotted
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if entry is invalid
  Assert(id != nullptr, MemoryException,WARNING,UNDEFINED);

  // Return and throw warning if memory store is empty
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  Assert(it != m_entry.end(), MemoryException,WARNING,NOT_FOUND);
  return (*it)->m_plot;
}

bool
Memory::getRestart(MemoryEntry* id)
//******************************************************************************
//  Return true if the variable is writted to restart file based on the ID
//! \param[in]  id  ID of the entry
//! \return True if the variable is written to restart file
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if entry is invalid
  Assert(id != nullptr, MemoryException,WARNING,UNDEFINED);

  // Return and throw warning if memory store is empty
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  Assert(it != m_entry.end(), MemoryException,WARNING,NOT_FOUND);
  return (*it)->m_restart;
}

const MemoryEntry*
Memory::getID(const string& name)
//******************************************************************************
//  Return the MemorySet key based on the variable name
//! \param[in]  name  Name of the variable
//! \return MemorySet key (used as ID)
//! \author J. Bakosi
//******************************************************************************
{
  Assert(name.size() != 0, MemoryException,FATAL,EMPTY_NAME);

  // Return and throw warning if memory store is empty
  Assert(m_name.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  // Find memory entry and return its value type
  auto it = m_name.find(name);
  Assert(it != m_name.end(), MemoryException,WARNING,NOT_FOUND);
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
  Assert(m_entry.size() != 0, MemoryException,WARNING,EMPTY_STORE);

  size_t bytes = 0;
  for (auto& e : m_entry) {
    bytes += sizeof(MemoryEntry) +
               e->m_number * VarComp[static_cast<int>(e->m_variable)] *
               SizeOf[static_cast<int>(e->m_value)];
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
  // Return and throw warning if entry is invalid
  Assert(id != nullptr, MemoryException,WARNING,UNDEFINED);

  // Get size of value type
  size_t size = SizeOf[static_cast<int>(id->m_value)];

  // Compute chunk size
  int i = id->m_number/m_nthreads;

  // Zero remaining portion
  memset(static_cast<char*>(id->m_ptr) + m_nthreads*i*size,
         0,
         (id->m_number%m_nthreads)*size);

  int myid;
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
