//******************************************************************************
/*!
  \file      src/Base/Memory.C
  \author    J. Bakosi
  \date      Sat 04 May 2013 07:45:54 AM MDT
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
#include <Paradigm.h>

using namespace Quinoa;

Memory::Memory(Paradigm* const paradigm) noexcept :
  m_nOMPthreads(paradigm->getOpenMP()->nthread())
//******************************************************************************
//  Constructor
//! \param[in]  paradigm Parallel programming object pointer
//! \author  J. Bakosi
//******************************************************************************
{
}

Memory::~Memory() noexcept
//******************************************************************************
//  Destructor
//! \details Free all allocated memory when leaving scope (just in case)
//! \author  J. Bakosi
//******************************************************************************
{
  freeAllEntries();
}

MemoryEntry*
Memory::newEntry(const size_t number,
                 const ValType value,
                 const VarType variable,
                 const string name,
                 const bool plot,
                 const bool restart)
//******************************************************************************
//  Allocate memory entry
//! \param[in]  number    Number of items to allocate
//! \param[in]  value     Chosen value type (BOOL_VAL, INT_VAL, REAL_VAL, etc.)
//! \param[in]  variable  Chosen variable type (SCALAR_VAR, VECTOR_VAR, etc.)
//! \param[in]  name      Given variable name
//! \param[in]  plot      True if variable can be plotted
//! \param[in]  restart   True if variable will be written to restart file
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the objects modified.
//! \return               ID for the allocated entry
//! \author     J. Bakosi
//******************************************************************************
{
  assert(number > 0);
  assert(name.size() > 0);

  // Compute total number of bytes to be allocated
  const size_t nbytes = number *
                        VarComp[static_cast<int>(variable)] *
                        SizeOf[static_cast<int>(value)];

  // Allocate memory
  void* ptr = static_cast<void*>(new (nothrow) char [nbytes]);
  if (ptr == nullptr)
    throw Exception(FATAL,
            "Cannot allocate memory for ptr in Memory::newEntry()");

  // Allocate memory entry metadata
  MemoryEntry* entry = new (nothrow) MemoryEntry(nbytes,
                                                 number,
                                                 value,
                                                 variable,
                                                 name,
                                                 plot,
                                                 restart,
                                                 ptr);
  if (entry == nullptr) {     // roll back changes and throw exception on error
    if (ptr) { delete [] static_cast<char*>(ptr); ptr = nullptr; }
    throw Exception(FATAL,
            "Cannot allocate memory for MemoryEntry in Memory::newEntry()");
  }

  // Store memory entry
  try {

    pair<MemorySet::iterator,bool> e = m_entry.insert(entry);
    if (!e.second) throw Exception(FATAL, "Cannot insert new memory entry");

  } // roll back changes and rethrow on error
    catch (exception&) {
      if (entry) { delete entry; entry = nullptr; }
      throw;
    }
    catch (...) {
      if (entry) { delete entry; entry = nullptr; }
      throw Exception(UNCAUGHT);
    }

  // Return key to caller
  return entry;
}

MemoryEntry*
Memory::newZeroEntry(const size_t number,
                     const ValType value,
                     const VarType variable,
                     const string name,
                     const bool plot,
                     const bool restart)
//******************************************************************************
//  Allocate and zero memory entry
//! \param[in]  number    Number of items to allocate
//! \param[in]  value     Chosen value type (BOOL_VAL, INT_VAL, REAL_VAL, etc.)
//! \param[in]  variable  Chosen variable type (SCALAR_VAR, VECTOR_VAR, etc.)
//! \param[in]  name      Given variable name
//! \param[in]  plot      True if variable can be plotted
//! \param[in]  restart   True if variable will be written to restart file
//! \return               ID for the allocated entry
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the objects modified.
//! \author     J. Bakosi
//******************************************************************************
{
  // Allocate new entry
  MemoryEntry* const entry =
    newEntry(number, value, variable, name, plot, restart);

  // Zero new entry
  zero(entry);

  // Return key to caller
  return entry;
}

void
Memory::freeEntry(MemoryEntry* id) noexcept
//******************************************************************************
//  Deallocate a memory entry
//! \param[in]  id  ID of the entry to be freed
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {

    assert(m_entry.size() != 0);

    // Return silently if entry is not yet allocated or already deallocated
    if (id == nullptr) return;

    // Find memory entry
    auto it = m_entry.find(id);
    assert(it != m_entry.end());

    // Deallocate memory entry pointed to by m_entry[id]
    // This also automatically calls MemoryEntry::~MemoryEntry(), which
    // deallocates the memory pointed to by MemoryEntry::m_ptr
    delete *it;

    // Remove MemoryEntry from MemorySet
    if (!m_entry.erase(id)) {
      throw Exception(WARNING, "Attempt to erase non-existent memory entry");
    }

    // Zero id, so the caller can tell the entry has been removed
    id = nullptr;

  } // emit only a warning on error
    catch (Exception& e) {
      e.echo("WARNING");
    }
    catch (exception& e) {
      cout << ">>> std::exception in Memory::freeEntry(): " << e.what() << endl;
    }
    catch (...) {
      cout << ">>> UNKNOWN EXCEPTION in Memory::freeEntry()" << endl;
    }
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
//!          exception here. Exception safety: no-throw guarantee: never throws
//!          exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  // Deallocate all memory entries
  // This also automatically calls MemoryEntry::~MemoryEntry(), which
  // deallocates the memory pointed to by MemoryEntry::m_ptr
  for (auto& e : m_entry) { delete e; }
  // Clear container
  m_entry.clear();
}

void
Memory::echoAllEntries(MemoryEntryField crit) const
//******************************************************************************
//  Echo all memory entries
//! \param[in]  crit    Sort items by memory entry field criteria
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if memory store is empty
  assert(m_entry.size() != 0);

  // Echo AllEntries-header
  cout << "* Dynamically allocated memory entries " << setfill(' ');
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
Memory::echo() const
//******************************************************************************
//  Echo unsorted memory entries
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

  for (auto* e : m_entry) cout << e->line();
}

void
Memory::echoByBytes() const
//******************************************************************************
//  Echo memory entries sorted by Bytes
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

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
Memory::echoByNumber() const
//******************************************************************************
//  Echo memory entries sorted by Number
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

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
Memory::echoByValue() const
//******************************************************************************
//  Echo memory entries sorted by Value
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

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
Memory::echoByVariable() const
//******************************************************************************
//  Echo memory entries sorted by Variable
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

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
Memory::echoByName() const
//******************************************************************************
//  Echo memory entries sorted by Name
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

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
Memory::echoByPlot() const
//******************************************************************************
//  Echo memory entries sorted by Plot
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

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
Memory::echoByRestart() const
//******************************************************************************
//  Echo memory entries sorted by Restart
//! \details    Exception safety: basic guarantee: if an exception is thrown,
//!             the stream remains in a valid state
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

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
Memory::getNumber(MemoryEntry* const id) const
//******************************************************************************
//  Return the number of items based on the ID
//! \param[in]  id  ID of the entry
//! \return     Number of items
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the memory entry container.
//! \author J. Bakosi
//******************************************************************************
{
  assert(id != nullptr);
  assert(m_entry.size() != 0);

  // Find memory entry and return its number of variables
  auto it = m_entry.find(id);
  assert(it != m_entry.end());
  return (*it)->m_number;
}

ValType
Memory::getValue(MemoryEntry* const id) const
//******************************************************************************
//  Return the value type based on the ID
//! \param[in]  id  ID of the entry
//! \return     The value type
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the memory entry container.
//! \author J. Bakosi
//******************************************************************************
{
  assert(id != nullptr);
  assert(m_entry.size() != 0);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  assert(it != m_entry.end());
  return (*it)->m_value;
}

VarType
Memory::getVariable(MemoryEntry* const id) const
//******************************************************************************
//  Return the variable type based on the ID
//! \param[in]  id  ID of the entry
//! \return     The variable type
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the memory entry container.
//! \author J. Bakosi
//******************************************************************************
{
  assert(id != nullptr);
  assert(m_entry.size() != 0);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  assert(it != m_entry.end());
  return (*it)->m_variable;
}

string
Memory::getName(MemoryEntry* const id) const
//******************************************************************************
//  Return the variable name based on the ID
//! \param[in]  id  ID of the entry
//! \return     The variable name
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the memory entry container.
//! \author J. Bakosi
//******************************************************************************
{
  assert(id != nullptr);
  assert(m_entry.size() != 0);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  assert(it != m_entry.end());
  return (*it)->m_name;
}

bool
Memory::getPlot(MemoryEntry* const id) const
//******************************************************************************
//  Return true if the variable can be plotted based on the ID
//! \param[in]  id  ID of the entry
//! \return     True if the variable can be plotted
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the memory entry container.
//! \author J. Bakosi
//******************************************************************************
{
  assert(id != nullptr);
  assert(m_entry.size() != 0);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  assert(it != m_entry.end());
  return (*it)->m_plot;
}

bool
Memory::getRestart(MemoryEntry* const id) const
//******************************************************************************
//  Return true if the variable is writted to restart file based on the ID
//! \param[in]  id  ID of the entry
//! \return     True if the variable is written to restart file
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the memory entry container.
//! \author J. Bakosi
//******************************************************************************
{
  assert(id != nullptr);
  assert(m_entry.size() != 0);

  // Find memory entry and return its value type
  auto it = m_entry.find(id);
  assert(it != m_entry.end());
  return (*it)->m_restart;
}

size_t
Memory::getBytes() const noexcept
//******************************************************************************
//  Return the number of allocated bytes
//! \details Return the number of bytes allocated in newEntry(). We account for
//!          the size of the MemoryEntry class instances and the allocated data
//!          pointed to by MemoryEntry::m_ptr. We do not account for the
//!          overhead of the STL container, therefore we will always
//!          underestimate the actual memory usage, though by only a very small
//!          fraction, i.e. <1e-4% for memory allocated in the range of MBytes.
//! \return  Number of allocated bytes in the memory store
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  assert(m_entry.size() != 0);

  size_t bytes = 0;
  for (auto& e : m_entry) {
    bytes += sizeof(MemoryEntry) +
               e->m_number * VarComp[static_cast<int>(e->m_variable)] *
               SizeOf[static_cast<int>(e->m_value)];
  }
  return bytes;
}

void
Memory::zero(MemoryEntry* const id) const noexcept
//******************************************************************************
//  Zero entry using multiple threads
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  assert(id != nullptr);

  // Get size of value type
  const size_t size = SizeOf[static_cast<int>(id->m_value)];

  // Compute chunk size
  const int i = id->m_number/m_nOMPthreads;

  // Zero remaining portion
  memset(static_cast<char*>(id->m_ptr) + m_nOMPthreads*i*size,
         0,
         (id->m_number%m_nOMPthreads)*size);

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
